clear all
close all
[p,t] = pear_mesh();

% Lineaire oplossing
% TODO: x0 = lineaire oplossing
x0 = 0;

% Volledige oplossing
% TODO: Zoek f en fp!
x = newton(x0,f,fp);