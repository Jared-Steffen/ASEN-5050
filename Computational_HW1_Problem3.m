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

% Givens
mu = 1;
r1 = [1,0,0];
r2 = [0,2,0];

% Calculate chord
c = r2 - r1;
