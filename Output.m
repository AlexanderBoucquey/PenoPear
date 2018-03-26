clear all;
close all;

p = load('meshp.txt');

Cu_0 = load('Cu_0.txt');
Cv_0 = load('Cv_0.txt');
Cu = load('Cu.txt');
Cv = load('Cv.txt');
% Plot lineaire oplossing oxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),Cu_0,xi,yi);
surf(xi,yi,zi);

% Plot lineaire oplossing dioxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),Cv_0,xi,yi);
surf(xi,yi,zi);

% Plot niet-lineaire oplossing oxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),Cu,xi,yi);
surf(xi,yi,zi);

% Plot niet-lineaire oplossing dioxide
figure
[xi,yi] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.1);
zi = griddata(p(:,1),p(:,2),Cv,xi,yi);
surf(xi,yi,zi);