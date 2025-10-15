clc; clear; close all

%% Problem Statement
%{
Write a different computer program “package” that will numerically
integrate the 2- body problem in Cartesian coordinates, i.e., solve the
equations of motion:

 ̇x = u
 ̇y = v
 ̇z = w
 ̇u = −μx/r3
 ̇v = −μy/r3
 ̇w = −μz/r3

where r = sqrt(x^2 + y^2 + z^2). The program should accept μ, initial
values of position and velocity, and a final time as input parameters.
For example, you can use the variable step RK4/5 integrator in Matlab, or
you can write your own. Make sure you have your error tolerances properly
set. Your program should:

(a) Produce data files suitable for plotting the trajectory
(x, y, z coordinates) and the velocity (u, v, w values) as a function of
time and against each other (i.e., as a 3-D plot showing the 
actual trajectory).

(b) Produce data files, suitable for plotting, containing the computed
values of the orbit elements, energy, and angular momentum at each point 
in time. For the 2-body problem, these should all be constant values.
Later we will use this program when we have additional forces, and we
will see that the orbit elements will then vary with time.

Use your programs from Question 1 to check your results. Specifically,
numerically integrate the orbit in Problem 1-i and compare with the
analytic solution you generated.
%}

% ICs
mu = 4e5; % km^3/s^2
r0 = [6 6 6] .* 1e3; % km
v0 = [-5 5 0]; % km/s
var = [r0';v0'];

% Determine orbital period for time vector
% Specific energy
epsilon = 0.5*norm(v0)^2-mu/norm(r0);
% Semi-major axis
a = -mu/(2*epsilon);
% Orbital Period
n = sqrt(mu/a^3);
T = (2*pi)/n;

% Time vector
t0 = 0;
t_step = 60;
tspan = t0:t_step:2*T;

% ode45 call
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t,state] = ode45(@(tspan,var) OrbitEOM(tspan,var,mu),tspan,var,options);

% Orbital Elements
[a,e,inc,Omega,w,tp,h_vec,epsilon] = orbital_elements(mu,state,t);

% Normalize h
h = zeros(height(state),1);
for j = 1:height(h_vec)
    h(j) = norm(h_vec(j,:));
end

% Plotting
ylabels_pos = ["$x$ [km]","$y$ [km]","$z$ [km]","$\dot{x}$ [km/s]",...
    "$\dot{y}$ [km/s]","$\dot{z}$ [km/s]"];

figure(1);
for i = 1:width(state)/2
    subplot(3,1,i)
    plot(t/3600,state(:,i),'LineWidth',2)
    grid on; grid minor
    xlabel('Time [hr]')
    ylabel(ylabels_pos(i),'Interpreter','latex')
end
sgtitle('Position Components for 2 Orbital Periods')

% exportgraphics(gcf, 'part_2_position_state_vars.png', 'Resolution', 300);

figure(2);
for i = 1:width(state)/2
    subplot(3,1,i)
    plot(t/3600,state(:,i+3),'LineWidth',2)
    grid on; grid minor
    xlabel('Time [hr]')
    ylabel(ylabels_pos(i+3),'Interpreter','latex')
end
sgtitle('Velocity Components for 2 Orbital Periods')

% exportgraphics(gcf, 'part_2_velocity_state_vars.png', 'Resolution', 300);

% Load in problem 1 results
problem1data = load("r_v_1i.mat");
r = problem1data.r;
v = problem1data.v;

figure(3);
plot3(state(:,1),state(:,2),state(:,3),'o')
hold on
plot3(r(1,:),r(2,:),r(3,:),'LineWidth',2)
grid on; grid minor
xlabel("$x$ [km]",'Interpreter','latex')
ylabel("$y$ [km]",'Interpreter','latex')
zlabel("$z$ [km]",'Interpreter','latex')
title('Position Space for 2 Orbital Periods')
legend('Integrated Results','Problem 1 Results','Location','best')

% exportgraphics(gcf, 'part_2_position_space.png', 'Resolution', 300);

figure(4);
plot3(state(:,4),state(:,5),state(:,6),'o')
hold on
plot3(v(1,:),v(2,:),v(3,:),'LineWidth',2)
grid on; grid minor
xlabel("$\dot{x}$ [km/s]",'Interpreter','latex')
ylabel("$\dot{y}$ [km/s]",'Interpreter','latex')
zlabel("$\dot{z}$ [km/s]",'Interpreter','latex')
title('Velocity Space for 2 Orbital Periods')
legend('Integrated Results','Problem 1 Results','Location','best')

% exportgraphics(gcf, 'part_2_velocity_space.png', 'Resolution', 300);

figure(5);
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

% exportgraphics(gcf, 'part_2_orbital_angles.png', 'Resolution', 300);

figure(6);
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

figure(7);
subplot(2,1,1)
plot(t/3600,h,'LineWidth',2)
xlabel('Time [hr]')
ylabel("Angular Momentum $h$ [$km^2/s$]",'Interpreter','latex')
grid on;grid minor
subplot(2,1,2)
plot(t/3600,epsilon,'LineWidth',2)
xlabel('Time [hr]')
ylabel("Energy $\varepsilon$",'Interpreter','latex')
grid on;grid minor
sgtitle('Integrals of Motion')

figure(8);
plot(t/3600,tp/3600,'LineWidth',2)
grid on; grid minor
xlabel('Time [hr]')
ylabel('Time of (Most Recent) Periapsis Passage [hr]')
title('Time of (Most Recent) Periapsis Passage')


% Output results
X = array2table(state,'VariableNames', {'x','y','z','u','v','w'});

vars = X.Properties.VariableNames;
for i = 1:numel(vars)
    if isnumeric(X.(vars{i}))
        X.(vars{i}) = round(X.(vars{i}), 5);
    end
end

writetable(X, 'state_data.csv');

oe_data = struct('a',a,'e',e,'i',rad2deg(inc),'w',rad2deg(w),'RAAN',rad2deg(Omega),'t_p',tp);
X = struct2table(oe_data);

vars = X.Properties.VariableNames;
for i = 1:numel(vars)
    if isnumeric(X.(vars{i}))
        X.(vars{i}) = round(X.(vars{i}), 5);
    end
end

writetable(X, 'oe_data.csv');

iom_data = struct('h',h,'Energy',epsilon);
X = struct2table(iom_data);

vars = X.Properties.VariableNames;
for i = 1:numel(vars)
    if isnumeric(X.(vars{i}))
        X.(vars{i}) = round(X.(vars{i}), 5);
    end
end

writetable(X, 'iom_data.csv');

%% Functions
function [var_dot] = OrbitEOM(~,var,mu)
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
    u_dot = (-mu*x)/r^3;
    v_dot = (-mu*y)/r^3;
    w_dot = (-mu*z)/r^3;

    % Final state derivative
    var_dot = [x_dot;y_dot;z_dot;u_dot;v_dot;w_dot];
end

function [a,e,inc,Omega,w,tp,h_vec,epsilon] = orbital_elements(mu,state,t)
    % Goal: Generate necessary orbital elements to describe an orbit
    
    % Extract r and v at each timestep
    r = state(:,1:3);
    v = state(:,4:6);

    % Unit vectors
    x_hat = [1 0 0]';
    y_hat = [0 1 0]';
    z_hat = [0 0 1]';

   % Initialize outputs
   h_vec = zeros(size(r));
   inc = zeros(height(state),1);
   a = zeros(height(state),1);
   e = zeros(height(state),1);
   w = zeros(height(state),1);
   Omega = zeros(height(state),1);
   epsilon = zeros(height(state),1);
   tp = zeros(height(state),1);

    for j = 1:height(r)
        % Angular momentum
        h_vec(j,:) = cross(r(j,:),v(j,:)); 
        h = norm(h_vec(j,:));
        h_hat = h_vec(j,:)./h;
    
        % Inclination
        inc(j) = acos(dot(h_hat,z_hat));
    
        % RAAN
        n_Omega_hat = cross(z_hat,h_hat)/norm(cross(z_hat,h_hat));
        Omega(j) = atan2(dot(n_Omega_hat,y_hat),dot(n_Omega_hat,x_hat));

        % Eccentricity
        e_vec = 1/mu.*cross(v(j,:),h_vec(j,:))-r(j,:)./norm(r(j,:));
        e(j) = norm(e_vec);
        e_hat = e_vec./e(j);
        e_hat_perp = cross(h_hat,e_hat);
    
        % Argument of periapsis
        n_Omega_hat_perp = cross(h_hat,n_Omega_hat);
        w(j) = atan2(dot(e_hat,n_Omega_hat_perp),dot(e_hat,n_Omega_hat));
    
        % Specific energy
        epsilon(j) = 0.5*norm(v(j,:))^2-mu/norm(r(j,:));
    
        % Semi-major axis
        a(j) = -mu/(2*epsilon(j));
    
        % Mean Motion
        n = sqrt(mu/a(j)^3);
    
        % Time of periapsis passage
        f = atan2(dot(r(j,:),e_hat_perp),dot(r(j,:),e_hat));
        E = 2*atan2(sqrt((1-e(j)))*sin(f/2),sqrt((1+e(j)))*cos(f/2));
        tp(j) = t(j) - (1/n)*mod((E-e(j)*sin(E)),2*pi);
    end
end