clc; clear; close all

%% Problem Statement
%{
Write a computer program “package” that will take an initial spacecraft 
position and velocity vector, r0, v0, at an initial time t0, and predict
the future/past position and velocity of the spacecraft at an arbitrary 
time t. To implement this you need to implement three distinct steps:

(a) Compute the orbital elements from an arbitrary position vector, 
velocity vector, and time.

(b) Solve Kepler's equation for an arbitrary time t to find the
corresponding value of true anomaly.

(c) Using the orbital elements, specify the position and velocity vectors
at the new value of true anomaly.

The program should take the gravitational parameter μ as an input.

i. Apply your program to the initial conditions:
μ = 4 × 10^5 km^3/s^2
r0 = [6 × 10^3, 6 × 10^3, 6 × 10^3] km
v0 = [−5, 5, 0] km/s

Propagate the orbit over two orbital periods using time intervals of
60 seconds. Plot the resulting trajectory in position and velocity space.

ii. Use the same value of μ as above. Assume an initial radius of periapsis
distance of 10,000 km, and orbit element angles i = 135◦, Ω = 45◦, and
ω = −90◦. For these given orbit elements, compute and plot the orbit in
x − y − z space over two orbit periods for: e = 0, 0.25, 0.5, 0.75, 0.9
%}

%% Part i)
% ICs
mu = 4e5; % km^3/s^2
r0 = [6 6 6] .* 1e3; % km
v0 = [-5 5 0]; % km/s
t0 = 0; % s

% Determine orbital elements
[a,p,e,i,Omega,w,T,n,tp,h,epsilon] = orbital_elements(mu,r0,v0,t0);

disp('Part 1i Results:')
disp("Semi-major axis: " + a + " km")
disp("Eccentricity: " + e)
disp("Inclination: " + rad2deg(i) + " degrees")
disp("RAAN: " + rad2deg(Omega) + " degrees")
disp("Argument of Periapsis: " + rad2deg(w) + " degrees")
disp("Time of Periapsis Passage (realtive to t0): " + tp + " s")
disp("Orbital Period: " + T + " s")
disp("Angular Momentum: " + h + "km^2/s")
disp("Energy: " + epsilon )

% Total simulation time
t_step = 60; % s
t = t0:t_step:2*T;

% Allowable error for Newton's iteration scheme for finding true anomaly
error = 1e-6;

% Solve Kepler's Equation
f = keplers_equation(e,t,tp,n,error);

% Solve 2BP
[r,v] = solution_2BP(mu,Omega,i,w,f,p,e);

% Save for comparison to problem 2
save("r_v_1i", "r", "v")

% Plot r and v in the position and velocity space
figure(1);
plot3(r(1,:),r(2,:),r(3,:),'LineWidth',2)
grid on; grid minor
xlabel('x distance [km]')
ylabel('y distance [km]')
zlabel('z distance [km]')
title('Position Space for 2 Orbital Periods')

% exportgraphics(gcf, 'part_1i_position_space.png', 'Resolution', 300);

figure(2);
plot3(v(1,:),v(2,:),v(3,:),'LineWidth',2)
grid on; grid minor
xlabel('x velocity component [km/s]')
ylabel('y velocity component [km/s]')
zlabel('z velocity component [km/s]')
title('Velocity Space for 2 Orbital Periods')

% exportgraphics(gcf, 'part_1i_velocity_space.png', 'Resolution', 300);

%% Part ii)

% Givens
rp = 10000; % km
i = deg2rad(135); % rad
Omega = deg2rad(45); % rad
w = deg2rad(-90); % rad
e = [0 0.25 0.5 0.75 0.99];

% Semi-major axis
a = rp./(1-e);

% Semi-latus rectum
p = a.*(1-e.^2);

% Time of periapsis passage and period
n = sqrt(mu./a.^3);
T = (2*pi)./n;
f0 = atan2(rp,0);
E0 = 2*atan2(sin(sqrt((1-e)/(1+e))*tan(f0/2)),cos(sqrt((1-e)/(1+e))*tan(f0/2)));
tp = t0 - (1./n).*(E0-e.*sin(E0));

% Total simulation time
t_step = 60; % s
for i = 1:length(T)    
    t = t0:t_step:2*T(i);
    f = keplers_equation(e,t,tp,n,error);
    [r,v] = solution_2BP(mu,Omega,i,w,f,p(i),e(i));
    figure(i+2);
    plot3(r(1,:),r(2,:),r(3,:),'LineWidth',2)
    grid on; grid minor
    xlabel('x distance [km]')
    ylabel('y distance [km]')
    zlabel('z distance [km]')
    title('Position Space for 2 Orbital Periods')
end


%% Functions
function [a,p,e,i,Omega,w,T,n,tp0,h,epsilon] = orbital_elements(mu,r0,v0,t0)
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
    tp0 = t0 - (1/n)*(E0-e*sin(E0));
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
                delta_E = -(M_star(j) - current_E+e(j)*sin(current_E))/...
                    (-1+e(j)*cos(current_E));
                current_E = current_E + delta_E;
                Euler_function = abs(M_star(j) - current_E+e(j)*sin(current_E));
            end
            f(i,j) = 2*atan2(sqrt(1+e(j))*sin(current_E/2),sqrt(1-e(j))*...
                cos(current_E/2));
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
    n_Omega_hat_perp = -cos(i)*sin(Omega)*x_hat + ...
        cos(i)*cos(Omega)*y_hat + sin(i)*z_hat;
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

