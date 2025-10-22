clc; clear; close all

%% part (c)(i)
% Problem: Plot inclinations as a function of orbit radius for circular Sun
% synchronous orbits

% Constants
mu = 4e5;
R_E = 6400;
J2 = 1e-3;
Omega_Sun = 2e-7;

% Maximum orbit radius
a_M = ((1.5*R_E^2*J2*sqrt(mu))/Omega_Sun)^(2/7); 

% Radius vector 
a = 6400:10:a_M;

% Inclination vector
i = acosd(-(a/a_M).^(7/2));

% Plot
figure();
plot(a,i,'LineWidth',2)
grid on; grid minor
xlabel('Orbit Radius [km]')
ylabel('Inclination [degrees]')
title('Possible Inclinations for Sun Synchronous Orbits')