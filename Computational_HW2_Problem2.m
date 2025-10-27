clc; clear; close all

%% Problem Statement
%{
Next, implement this J2 problem in your numerical integration code to check
and verify your analytical solutions, and vice-versa.

(a) The acceleration of the J2 effect acting on an orbit is found from the
full force potential

uJ2 = (−3*μ*R0^2*J2)/(2*r^7) * [[x^2 + y^2 − 4z^2]x,
[x^2 + y^2 − 4z^2]y, [3(x^2 + y^2) − 2z^2]z]

Verify this relationship by taking the gradient of the initial RJ2 potential.
Implement this perturbation acceleration in your numerical integration code.
To verify your numerical integration, show that the following integrals are
constant in your simulations.

Angular momentum about the z-axis: hz = z_hat*(r × v).
Energy: E = 0.5*v^2-μ/r-μ/(2r^3)*R0^2*J2*[1-3*(z^2/r^2)]

(b) For a given set of initial orbit elements of a = 7000, e = 0.1, i = 45◦
and ω = Ω = σ = 0, compare your numerical integration to the predicted orbit
elements from the averaging theory over at least 5 orbits. Do this by
plotting the orbit elements from the numerical solution and the averaged
solution as functions of time on the same graph.

c) Note that the average rates of the argument of periapsis and longitude of
the ascending node will diverge as time goes on. This is due to the short
period terms – specifically, the true "average” semi-major axis, eccentricity
and inclination are different than the initial conditions that you used.

Adjust the initial conditions of the averaged solution to equal the observed
average values of your numerically computed a, e and i. Does the agreement
between the average rates of argument of periapsis and longitude of the
ascending node agree more precisely?

%}

% Constants
mu = 4e5; % km^3/s^2
RE = 6400; % km
J2 = 1e-3;

%% part (a)

% ICs
r0 = [6 6 6] .* 1e3; % km
v0 = [-5 5 0]; % km/s
var = [r0';v0'];

% Determine orbital period for time vector
% Specific energy
epsilon = 0.5*norm(v0)^2-mu/norm(r0)...
    -mu/(2*norm(r0)^3)*RE^2*J2*(1-3*r0(3)^2/norm(r0)^2);
% Semi-major axis
a = -mu/(2*epsilon);
% Orbital Period
n = sqrt(mu/a^3);
T = (2*pi)/n;

% Time vector
t0 = 0;
t_step = 60;
tspan = t0:t_step:T;

% ode45 call
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode45(@(tspan,var) OrbitEOM_J2(tspan,var,mu,RE,J2),tspan,var,options);

% Orbital Elements
[~,~,~,~,~,~,hz,epsilon] = orbital_elements_J2(mu,RE,J2,state,t);

figure(1);
subplot(2,1,1)
plot(t/3600,hz,'LineWidth',2)
xlabel('Time [hr]')
ylabel("Angular Momentum about the z-axis $h_z$ [$km^2/s$]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t/3600,epsilon,'LineWidth',2)
xlabel('Time [hr]')
ylabel("Energy $\varepsilon$",'Interpreter','latex')
grid on;grid minor
sgtitle('Integrals of Motion')

%% part (b)

% Givens
a = 7000; % km
e = 0.1;
i = deg2rad(45); % rad
w = 0; % rad
Omega = 0; % rad
tp = 0; % rad

% Finding intial r and v (assume at periapsis)
rp = a*(1-e);
vp = sqrt(mu/a*((1+e)/(1-e)));

% Unit vectors
x_hat = [1 0 0]';
y_hat = [0 1 0]';
z_hat = [0 0 1]';
n_Omega_hat = cos(Omega)*x_hat + sin(Omega)*y_hat;
n_Omega_hat_perp = -cos(i)*sin(Omega)*x_hat + ...
    cos(i)*cos(Omega)*y_hat + sin(i)*z_hat;
e_hat = cos(w)*n_Omega_hat + sin(w)*n_Omega_hat_perp;
e_hat_perp = -sin(w)*n_Omega_hat + cos(w)*n_Omega_hat_perp;

r0 = rp*e_hat;
v0 = vp*e_hat_perp;
var = [r0;v0];

% Orbital Period
n = sqrt(mu/a^3);
T = (2*pi)/n;

% Time vector
t_step = 60;
tspan = tp:t_step:10*T;

% ode45 call
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state2] = ode45(@(tspan,var) OrbitEOM_J2(tspan,var,mu,RE,J2),tspan,var,options);

% Orbital Elements
[a,e,inc,Omega,w,sigma,hz,epsilon] = orbital_elements_J2(mu,RE,J2,state2,t);

figure(2);
subplot(2,1,1)
plot(t/3600,a,'LineWidth',2)
xlabel('Time [hr]')
ylabel("Semi-major axis $a$ [km]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t/3600,e,'LineWidth',2)
xlabel('Time [hr]')
ylabel("Eccetricity $e$",'Interpreter','latex')
grid on;grid minor
sgtitle('Orbital Elements: Geometry')

figure(3);
subplot(3,1,1)
plot(t/3600,rad2deg(Omega),'LineWidth',2)
xlabel('Time [hr]')
ylabel("$\Omega$ [degrees]",'Interpreter','latex')
grid on;grid minor
subplot(3,1,2)
plot(t/3600,rad2deg(w),'LineWidth',2)
ylabel("$\omega$ [degrees]",'Interpreter','latex')
grid on;grid minor
subplot(3,1,3)
plot(t/3600, rad2deg(inc),'LineWidth',2)
ylabel("$i$ [degrees]",'Interpreter','latex')
grid on;grid minor
sgtitle('Orbital Elements: Angles')

figure(4);
plot(t/3600,sigma,'LineWidth',2)
grid on; grid minor
xlabel('Time [hr]')
ylabel('$\sigma$ [degrees]','Interpreter','latex')
title('Orbital Elements: $\sigma = -nt_p$','Interpreter','latex')



%% part (c)

%% Functions
function [var_dot] = OrbitEOM_J2(~,var,mu,RE,J2)
    % Goal: Output ODEs for ode45

    % Extract state variables
    x = var(1);
    y = var(2);
    z = var(3);
    u = var(4);
    v = var(5);
    w = var(6);

    % Calculate radius
    r = sqrt(x^2 + y^2 + z^2);

    % J2 Acceleration Coefficient
    uJ2_coeff = (-1.5*mu*RE^2*J2)/r^7;

    % Assign t.r.o.c variables
    x_dot = u;
    y_dot = v;
    z_dot = w;
    u_dot = (-mu*x)/r^3 + uJ2_coeff*(x^2+y^2-4*z^2)*x;
    v_dot = (-mu*y)/r^3 + uJ2_coeff*(x^2+y^2-4*z^2)*y;
    w_dot = (-mu*z)/r^3 + uJ2_coeff*(3*(x^2+y^2)-2*z^2)*z;

    % Final state derivative
    var_dot = [x_dot;y_dot;z_dot;u_dot;v_dot;w_dot];
end

function [a,e,inc,Omega,w,sigma,hz_vec,epsilon] = orbital_elements_J2(mu,RE,J2,state,t)
    % Goal: Generate necessary orbital elements to describe an orbit
    
    % Extract r and v at each timestep
    r = state(:,1:3);
    v = state(:,4:6);

    % Unit vectors
    x_hat = [1 0 0]';
    y_hat = [0 1 0]';
    z_hat = [0 0 1]';

   % Initialize outputs
   hz_vec = zeros(height(state),1);
   inc = zeros(height(state),1);
   a = zeros(height(state),1);
   e = zeros(height(state),1);
   w = zeros(height(state),1);
   Omega = zeros(height(state),1);
   epsilon = zeros(height(state),1);
   sigma = zeros(height(state),1);

    for j = 1:height(r)
        % Angular momentum
        h_vec = cross(r(j,:),v(j,:)); 
        h = norm(h_vec);
        h_hat = h_vec./h;
        hz_vec(j,:) = dot(z_hat, h_vec);
    
        % Inclination
        inc(j) = acos(dot(h_hat,z_hat));
    
        % RAAN
        n_Omega_hat = cross(z_hat,h_hat)/norm(cross(z_hat,h_hat));
        Omega(j) = atan2(dot(n_Omega_hat,y_hat),dot(n_Omega_hat,x_hat));

        % Eccentricity
        e_vec = 1/mu.*cross(v(j,:),h_vec)-r(j,:)./norm(r(j,:));
        e(j) = norm(e_vec);
        e_hat = e_vec./e(j);
        e_hat_perp = cross(h_hat,e_hat);
    
        % Argument of periapsis
        n_Omega_hat_perp = cross(h_hat,n_Omega_hat);
        w(j) = atan2(dot(e_hat,n_Omega_hat_perp),dot(e_hat,n_Omega_hat));
    
        % Specific energy
        epsilon(j) = 0.5*norm(v(j,:))^2-mu/norm(r(j,:))...
            -mu/(2*norm(r(j,:))^3)*RE^2*J2*(1-3*r(j,3)^2/norm(r(j,:))^2);
    
        % Semi-major axis
        a(j) = h^2/(mu*(1-e(j)^2));
        
        % Mean Motion
        n = sqrt(mu/a(j)^3);
    
        % Time of periapsis passage
        f = atan2(dot(r(j,:),e_hat_perp),dot(r(j,:),e_hat));
        E = 2*atan2(sqrt((1-e(j)))*sin(f/2),sqrt((1+e(j)))*cos(f/2));
        tp = t(j) - (1/n)*(E - e(j)*sin(E));
        sigma(j) = -n*tp;
    end
end
