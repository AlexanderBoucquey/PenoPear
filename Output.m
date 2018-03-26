clear all;
close all;

Cu_0 = load('Cu_0.txt');
Cv_0 = load('Cv_0.txt');
Cu = load('Cu.txt');
Cv = load('Cv.txt');
% Plot lineaire oplossing oxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),C_u0,xi,yi);
surf(xi,yi,zi);

% Plot lineaire oplossing dioxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),C_v0,xi,yi);
surf(xi,yi,zi);

figure
scatter(rp(:,1),rp(:,2));