clc; clear; close all

%% Problem Statement
%{
Next, implement the constant acceleration problem in your numerical
integration code to check and verify your analytical solutions, and vice-versa.

(a) Numerically simulate an initially circular trajectory for a value of
μ = 1, a = 1, i = 0 and g = 0.001, 0.01, 0.1, 1.0. Use the following acceleration
as a perturbation to your 2-body numerical integration code.

u = [g; 0; 0]

To verify your numerical integration, show that the following integrals
are constant in your simulations.

Angular momentum about the x-axis: hx = x_hat*(r × v).
Energy: E = 0.5*v^2-μ/r-gx

(b) Plot the trajectory and the orbit elements of a, e  ̃ω for the different
cases of acceleration and compare with your analytical solution for the
averaged case. Discuss the comparison.

(c) Is it possible to improve the agreement by making the short-period shift
used above for the J2 problem? What are the limits on the acceleration magnitude
g for the averaged solution to reasonably model the motion?
%}

% Givens/ICs
a = 1; 
mu = 1;
i = 0;
g = [0.001,0.01,0.1,1.0];

% Assuming starting along z-axis
r = a;
r0 = [0 0 r];
v0 = [0 sqrt(mu/r) 0];
var = [r0',v0'];

% Orbital Period
n = sqrt(mu/a^3);
T = (2*pi)/n;

% Time vector
t0 = 0; % s
t_step = 0.1; % s
tspan = t0:t_step:2*T;

% ode45 call
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t1,state1] = ode45(@(tspan,var) OrbitEOM_constF(tspan,var,mu,g(1)),tspan,var,options);
[t2,state2] = ode45(@(tspan,var) OrbitEOM_constF(tspan,var,mu,g(2)),tspan,var,options);
[t3,state3] = ode45(@(tspan,var) OrbitEOM_constF(tspan,var,mu,g(3)),tspan,var,options);
[t4,state4] = ode45(@(tspan,var) OrbitEOM_constF(tspan,var,mu,g(4)),tspan,var,options);

% Orbital Elements
[a1,e1,w_tilde1,hx1,epsilon1] = orbital_elements_constF(mu,g(1),state1);
[a2,e2,w_tilde2,hx2,epsilon2] = orbital_elements_constF(mu,g(2),state2);
[a3,e3,w_tilde3,hx3,epsilon3] = orbital_elements_constF(mu,g(3),state3);
[a4,e4,w_tilde4,hx4,epsilon4] = orbital_elements_constF(mu,g(4),state4);

%% Plots
figure(1);
subplot(2,1,1)
plot(t1,hx1,'LineWidth',2)
xlabel('Time [s]')
ylabel("Angular Momentum about the x-axis $h_x$ [$km^2/s$]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t1,epsilon1,'LineWidth',2)
xlabel('Time [s]')
ylabel("Energy $\varepsilon$",'Interpreter','latex')
grid on;grid minor
sgtitle('Integrals of Motion for g = 0.001')

figure(2);
subplot(2,1,1)
plot(t2,hx2,'LineWidth',2)
xlabel('Time [s]')
ylabel("Angular Momentum about the x-axis $h_x$ [$km^2/s$]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t2,epsilon2,'LineWidth',2)
xlabel('Time [s]')
ylabel("Energy $\varepsilon$",'Interpreter','latex')
grid on;grid minor
sgtitle('Integrals of Motion for g = 0.01')

figure(3);
subplot(2,1,1)
plot(t3,hx3,'LineWidth',2)
xlabel('Time [s]')
ylabel("Angular Momentum about the x-axis $h_x$ [$km^2/s$]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t3,epsilon3,'LineWidth',2)
xlabel('Time [s]')
ylabel("Energy $\varepsilon$",'Interpreter','latex')
grid on;grid minor
sgtitle('Integrals of Motion for g = 0.1')

figure(4);
subplot(2,1,1)
plot(t4,hx4,'LineWidth',2)
xlabel('Time [s]')
ylabel("Angular Momentum about the x-axis $h_x$ [$km^2/s$]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t4,epsilon4,'LineWidth',2)
xlabel('Time [s]')
ylabel("Energy $\varepsilon$",'Interpreter','latex')
grid on;grid minor
sgtitle('Integrals of Motion for g = 1.0')


%% Functions
function [var_dot] = OrbitEOM_constF(~,var,mu,g)
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

    % Assign t.r.o.c variables
    x_dot = u;
    y_dot = v;
    z_dot = w;
    u_dot = (-mu*x)/r^3 + g;
    v_dot = (-mu*y)/r^3;
    w_dot = (-mu*z)/r^3;

    % Final state derivative
    var_dot = [x_dot;y_dot;z_dot;u_dot;v_dot;w_dot];
end

function [a,e,w_tilde,hx_vec,epsilon] = orbital_elements_constF(mu,g,state)
    % Goal: Generate necessary orbital elements to describe an orbit
    
    % Extract r and v at each timestep
    r = state(:,1:3);
    v = state(:,4:6);
    x = r(:,1);

    % Unit vectors
    x_hat = [1 0 0]';
    y_hat = [0 1 0]';
    z_hat = [0 0 1]';

   % Initialize outputs
   hx_vec = zeros(height(state),1);
   a = zeros(height(state),1);
   e = zeros(height(state),1);
   epsilon = zeros(height(state),1);
   w_tilde = zeros(height(state),1);

    for j = 1:height(r)
        % Angular momentum
        h_vec = cross(r(j,:),v(j,:)); 
        h = norm(h_vec);
        h_hat = h_vec./h;
        hx_vec(j,:) = dot(x_hat, h_vec);

        % RAAN
        n_Omega_hat = cross(z_hat,h_hat)/norm(cross(z_hat,h_hat));
        Omega = atan2(dot(n_Omega_hat,y_hat),dot(n_Omega_hat,x_hat));

        % Eccentricity
        e_vec = 1/mu.*cross(v(j,:),h_vec)-r(j,:)./norm(r(j,:));
        e(j) = norm(e_vec);
        e_hat = e_vec./e(j);

        % Argument of periapsis
        n_Omega_hat_perp = cross(h_hat,n_Omega_hat);
        w = atan2(dot(e_hat,n_Omega_hat_perp),dot(e_hat,n_Omega_hat));

        % Longitude of Peripasis
        w_tilde(j) = w + Omega;
    
        % Specific energy
        epsilon(j) = 0.5*norm(v(j,:))^2-mu/norm(r(j,:))-g*x(j);
    
        % Semi-major axis
        a(j) = h^2/(mu*(1-e(j)^2));
       
    end
end