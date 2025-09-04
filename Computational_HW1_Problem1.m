clc; clear; close all

%% Problem Statement
%{
Write a computer program “package” that will take an initial spacecraft 
position and velocity vector, r_o, v_o, at an initial time t_o, and predict
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
r_o = [6 × 10^3, 6 × 10^3, 6 × 10^3] km
v_o = [−5, 5, 0] km/s

Propagate the orbit over two orbital periods using time intervals of
60 seconds.Plot the resulting trajectory in position and velocity space.

ii. Use the same value of μ as above. Assume an initial radius of periapsis
distance of 10,000 km, and orbit element angles i = 135◦, Ω = 45◦, and
ω = −90◦. For these given orbit elements, compute and plot the orbit in
x − y − z space over two orbit periods for: e = 0, 0.25, 0.5, 0.75, 0.9
%}

% ICs
mu = 4e5; % km^3/s^2
r_o = [6 6 6] .* 1e3; % km
v_o = [-5 5 0]; % km/s


function[x,y,z] = orbital_elements(mu,r_o,v_o)

end

