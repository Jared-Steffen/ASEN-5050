clc; clear; close all

%% Problem Statement
%{
Write a third computer program “package” that will solve Lambert’s problem: Given
two position vectors and a specified interval of time, compute the necessary initial
velocity. Your program should:

(a) Determine if the requested transfer is a hyperbolic or elliptic orbit.

(b) Compute both elliptic transfers (θ < π, θ > π), if they exist.

Use your programs from Questions 1 and 2 to check your results. Then carry out the
following computations, all using μ = 1 and a given set of position vectors r1 = [1, 0, 0]
and r2 = [0, 2, 0].

(a) What is the minimum energy transfer ellipse between these two points? What
times of transfer correspond to this transfer ellipse?

(b) What is the minimum eccentricity transfer ellipse between these two points? What
times of transfer correspond to this transfer ellipse?

(c) What is the parabolic transfer time between these two points? Choose transfer
times equal to twice of the parabolic transfer times and solve Lambert’s problem
between these two points.
%}

%% Testing Algorithm w/ Problem 1 and 2 Code
% Givens/ICs
mu1 = 4e5; % km^3/s^2
r1_vec = [5,5,5].*1e3; % km
r2_vec = [-5,-6,5].*1e3; % km
delta_t = 1800; % s

[v1,~,~,~,~] = LambertsTheorem(mu1,delta_t,r1_vec,r2_vec);

% Use functions from problem 1 to verify both transfer ellipses
t0 = 0; % s
[a1,p1,e1,i1,Omega1,w1,T1,n1,tp1] = orbital_elements(mu1,r1_vec,v1(1,:),t0);
[a2,p2,e2,i2,Omega2,w2,T2,n2,tp2] = orbital_elements(mu1,r1_vec,v1(2,:),t0);
t_step = 60; % s
t1 = t0:t_step:T1;
t2 = t0:t_step:T2;

error = 1e-6;
f1 = keplers_equation(e1,t1,tp1,n1,error);
f2 = keplers_equation(e2,t2,tp2,n2,error);

[r1_1,~] = solution_2BP(mu1,Omega1,i1,w1,f1,p1,e1);
[r1_2,~] = solution_2BP(mu1,Omega2,i2,w2,f2,p2,e2);

% Plot r and v in the position and velocity space
figure(1);
plot3(r1_1(1,:),r1_1(2,:),r1_1(3,:),'LineWidth',2)
hold on
plot3(r1_2(1,:),r1_2(2,:),r1_2(3,:),'LineWidth',2)
plot3(r1_vec(1),r1_vec(2),r1_vec(3),'.','MarkerSize',50)
plot3(r2_vec(1),r2_vec(2),r2_vec(3),'.','MarkerSize',50)
grid on; grid minor
xlabel("$x$ [km]",'Interpreter','latex')
ylabel("$y$ [km]",'Interpreter','latex')
zlabel("$z$ [km]",'Interpreter','latex')
title("Transfer Orbits for Lambert's Problem")
legend('Transfer Orbit 1','Transfer Orbit 2','r1','r2')

% Use functions from problem 2 to verify both transfer ellipses
var1 = [r1_vec';v1(1,:)'];
var2 = [r1_vec';v1(2,:)'];
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t1,state1] = ode45(@(t1,var1) OrbitEOM(t1,var1,mu1),t1,var1,options);
[t2,state2] = ode45(@(t2,var2) OrbitEOM(t1,var2,mu1),t2,var2,options);

figure(2);
plot3(state1(:,1),state1(:,2),state1(:,3),'LineWidth',2)
hold on
plot3(state2(:,1),state2(:,2),state2(:,3),'LineWidth',2)
plot3(r1_vec(1),r1_vec(2),r1_vec(3),'.','MarkerSize',50)
plot3(r2_vec(1),r2_vec(2),r2_vec(3),'.','MarkerSize',50)
grid on; grid minor
xlabel("$x$ [km]",'Interpreter','latex')
ylabel("$y$ [km]",'Interpreter','latex')
zlabel("$z$ [km]",'Interpreter','latex')
title("Transfer Orbits for Lambert's Problem")
legend('Transfer Orbit 1','Transfer Orbit 2','r1','r2')

%% Calculations for minimum (energy, eccentricity) transfer ellipses
% Givens/ICs
mu2 = 1; % km^3/s^2
r1_2_vec = [1,0,0]; % km
r2_2_vec = [0,2,0]; % km
delta_t1 = 3.7712; % s
delta_t2 = 4.2164; % s

[v1_2a,a_m,t_m,e_m,t_pm] = LambertsTheorem(mu2,delta_t1,r1_2_vec,r2_2_vec);
[v1_2b,a_m,t_m,e_m,t_pm] = LambertsTheorem(mu2,delta_t2,r1_2_vec,r2_2_vec);

% Minimum Eccentricity TOF
syms w
eqn1 = norm(r1_2_vec)*(1+e_m*cos(0-w)) == norm(r2_2_vec)*(1+e_m*cos((pi/2)-w));

omega = double(vpasolve(eqn1, w));
f1 =  0 - omega;
f2 = pi/2 - omega;

E1 = 2*atan2(sqrt((1-e_m))*sin(f1/2),sqrt((1+e_m))*cos(f1/2));
E2 = 2*atan2(sqrt((1-e_m))*sin(f2/2),sqrt((1+e_m))*cos(f2/2));

p = norm(r1_2_vec)*(1+e_m*cos(0-omega));
a = p/(1-e_m^2);
T = (2*pi)/sqrt(mu2)*a^(1.5);

ToF1 = (a^(1.5)*(E2-E1-e_m*(sin(E2)-sin(E1))))/sqrt(mu2);
ToF2 = T - ToF1;

% 2*Parabolic ToF Plot (using problem 1 and 2 code)
t0 = 0; % s
[a1,p1,e1,i1,Omega1,w1,T1,n1,tp1] = orbital_elements(mu2,r1_2_vec,v1_2a(1,:),t0);
[a2,p2,e2,i2,Omega2,w2,T2,n2,tp2] = orbital_elements(mu2,r1_2_vec,v1_2b(1,:),t0);
t_step = 0.01; % s
t1 = t0:t_step:T1;
t2 = t0:t_step:T2;

error = 1e-6;
f1 = keplers_equation(e1,t1,tp1,n1,error);
f2 = keplers_equation(e2,t2,tp2,n2,error);

[r1_1,~] = solution_2BP(mu2,Omega1,i1,w1,f1,p1,e1);
[r1_2,~] = solution_2BP(mu2,Omega2,i2,w2,f2,p2,e2);

var1 = [r1_2_vec';v1_2a(1,:)'];
var2 = [r1_2_vec';v1_2b(1,:)'];
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t1,state1] = ode45(@(t1,var1) OrbitEOM(t1,var1,mu2),t1,var1,options);
[t2,state2] = ode45(@(t2,var2) OrbitEOM(t1,var2,mu2),t2,var2,options);

figure(3);
plot3(state1(:,1),state1(:,2),state1(:,3),'LineWidth',2)
hold on
plot3(state2(:,1),state2(:,2),state2(:,3),'LineWidth',2)
plot3(r1_2_vec(1),r1_2_vec(2),r1_2_vec(3),'.','MarkerSize',50)
plot3(r2_2_vec(1),r2_2_vec(2),r2_2_vec(3),'.','MarkerSize',50)
grid on; grid minor
xlabel("$x$ [km]",'Interpreter','latex')
ylabel("$y$ [km]",'Interpreter','latex')
zlabel("$z$ [km]",'Interpreter','latex')
title("Transfer Orbits for Lambert's Problem when $\Delta t = 2t_{parabolic}$",'Interpreter','latex')
legend('Transfer Orbit 1','Transfer Orbit 2','r1','r2')

%% Functions (some are taken from problems 1 and 2)
function [v1,a_m,t_m,e_m,tp_m] = LambertsTheorem(mu,delta_t,r1_vec,r2_vec)

% Calculate magnitude of distances and angle between
r1 = norm(r1_vec);
r1_hat = r1_vec./r1;
r2 = norm(r2_vec);
theta1 = acos(dot(r1_vec,r2_vec)/(r2*r1));
theta2 = 2*pi - theta1;
theta = [theta1,theta2];

% Calculate chord
c_vec = r2_vec - r1_vec;
c = norm(c_vec);
c_hat = c_vec./c;

% Minimum energy ellipse and TOFs
a_m = 0.25*(r1+r2+c);
s = 0.5*(r1+r2+c);
beta_m = 2*asin(sqrt((s-c)/s));
t_m = ((pi-beta_m+sin(beta_m))*sqrt(s^3/8))/sqrt(mu);
n_m = sqrt(mu/a_m^3);
T_m = (2*pi)/n_m;

% Determine both travels times per ellipse
t_F2 = T_m-delta_t;
t_F = [delta_t,t_F2];

% Minimum eccentricity ellipse
e_m = abs((r2-r1))/c;

% % Solve for semi-major axis, alpha, beta
syms a alpha0 beta0 real

v1 = [];
tp_m = [];
for i = 1:length(theta)
    if sin(theta(i)) > 0
        tp = (sqrt(2)/(3*sqrt(mu)))*(s^(3/2)-(s-c)^(1.5));
        beta2 = beta0;
        tp_m = [tp_m;tp];
    else
        tp = (sqrt(2)/(3*sqrt(mu)))*(s^(3/2)+(s-c)^(1.5));
        beta2 = -beta0;
        tp_m = [tp_m;tp];
    end

    if t_F(i) <= t_m
        alpha2 = alpha0;
    else
        alpha2 = 2*pi - alpha0;
    end
    eqn1 = a^1.5*(alpha2-beta2-(sin(alpha2)-sin(beta2))) == sqrt(mu)*delta_t;
    eqn2 = 1/(4*a)*(r1+r2+c) == sin(alpha2/2)^2;
    eqn3 = 1/(4*a)*(r1+r2-c) == sin(beta2/2)^2;

    soln = vpasolve([eqn1,eqn2,eqn3],[a,alpha0,beta0]);
    a2 = double(soln.a);
    a2_val = real(a2^2);
    alpha = double(soln.alpha0);
    beta = double(soln.beta0);
    if a2_val < 0
        disp("This transfer is hyperbolic")
        v1 = [v1;NaN,NaN,NaN];
    else
        disp("This transfer is elliptic")
        A = sqrt(mu/(4*a2))*cot(alpha/2);
        B = sqrt(mu/(4*a2))*cot(beta/2);
        v1 = [v1;(B+A).*c_hat + (B-A).*r1_hat];
    end
end
end



function [a,p,e,i,Omega,w,T,n,tp] = orbital_elements(mu,r0,v0,t0)
    % Goal: Generate necessary orbital elements to describe an orbit

    % Unit vectors
    x_hat = [1 0 0]';
    y_hat = [0 1 0]';
    z_hat = [0 0 1]';

    % Angular momentum
    h_vec = cross(r0,v0); 
    h = norm(h_vec);
    h_hat = h_vec./h;
    
    % Inclination
    i = acos(dot(h_hat,z_hat));

    % Semi-latus rectum
    p = h^2/mu;

    % RAAN
    n_Omega_hat = cross(z_hat,h_hat)/norm(cross(z_hat,h_hat));
    Omega = atan2(dot(n_Omega_hat,y_hat),dot(n_Omega_hat,x_hat));

    % Eccentricity
    e_vec = 1/mu*cross(v0,h_vec)-r0/norm(r0);
    e = norm(e_vec);
    e_hat = e_vec./e;
    e_hat_perp = cross(h_hat,e_hat);

    % Argument of periapsis
    n_Omega_hat_perp = cross(h_hat,n_Omega_hat);
    w = atan2(dot(e_hat,n_Omega_hat_perp),dot(e_hat,n_Omega_hat));

    % Specific energy
    epsilon = 0.5*norm(v0)^2-mu/norm(r0);

    % Semi-major axis
    a = -mu/(2*epsilon);

    % Orbital Period
    n = sqrt(mu/a^3);
    T = (2*pi)/n;

    % Time of periapsis passage
    f0 = atan2(dot(r0,e_hat_perp),dot(r0,e_hat));
    E0 = 2*atan2(sqrt((1-e))*sin(f0/2),sqrt((1+e))*cos(f0/2));
    tp = t0 - (1/n)*(E0-e*sin(E0));
end

function f = keplers_equation(e,t,tp,n,error)
    % Goal: Generate true anomaly at each time step in the simulation

    % Initialize true anomaly
    f = zeros(length(t),length(e));

    % Iteration scheme for each time step
    for j = 1:length(e)
        for i = 1:length(t)
            M_star = n*(t(i)-tp(j));
            current_E = M_star(j);
            Euler_function = abs(M_star(j) - current_E+e(j)*sin(current_E));
            while Euler_function > error
                delta_E = -(M_star(j) - current_E+e(j)*sin(current_E))/(-1+e(j)*cos(current_E));
                current_E = current_E + delta_E;
                Euler_function = abs(M_star(j) - current_E+e(j)*sin(current_E));
            end
            f(i,j) = 2*atan2(sqrt(1+e(j))*sin(current_E/2),sqrt(1-e(j))*cos(current_E/2));
        end
    end
end


function [r,v] = solution_2BP(mu,Omega,i,w,f,p,e)
    % Goal: Determine r and v at every time step of the simulation

    % Unit vectors
    x_hat = [1 0 0]';
    y_hat = [0 1 0]';
    z_hat = [0 0 1]';
    n_Omega_hat = cos(Omega)*x_hat + sin(Omega)*y_hat;
    n_Omega_hat_perp = -cos(i)*sin(Omega)*x_hat + cos(i)*cos(Omega)*y_hat + sin(i)*z_hat;
    e_hat = cos(w)*n_Omega_hat + sin(w)*n_Omega_hat_perp;
    e_hat_perp = -sin(w)*n_Omega_hat + cos(w)*n_Omega_hat_perp;

    % Initialize r and v vectors
    r = zeros(3,length(f));
    v = zeros(3,length(f));

    % r and v for each true anomaly
    for i = 1:length(f)
        r(:,i) = p/(1+e*cos(f(i))) * (cos(f(i))*e_hat + sin(f(i))*e_hat_perp);
        v(:,i) = sqrt(mu/p) * (-sin(f(i))*e_hat + (e+cos(f(i)))*e_hat_perp);
    end

end

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