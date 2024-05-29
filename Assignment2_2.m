clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% GEOID
file = 'JGDWN_CER18D_GEOID_0018.IMG';
res = 1;

fid = fopen(file,'r','ieee-be');
data = fread(fid,[360 180])';
fclose(fid); 

latLim = [-89.5 89.5 res];
lonLim = [-179.5 179.5 res];
latGrid = latLim(1):latLim(3):latLim(2);
lonGrid = lonLim(1):lonLim(3):lonLim(2);

data = [data(:, 181:end), data(:, 1:180)];

a = 482e3;
f = 0.074896;
e2 = 2*f - f^2;
[~, latMesh] = meshgrid(lonGrid, latGrid);
phi = deg2rad(latMesh);
great_normal = (a ./ sqrt(1 - e2 * sin(phi).^2));
x_ellipsoid = great_normal .* cos(phi);
y_ellipsoid = great_normal .* sin(phi) * (1 - e2);
ellipsoid = sqrt(x_ellipsoid.^2 + y_ellipsoid.^2);

geoid = ellipsoid + data;

figure;
imagesc(lonGrid, latGrid, geoid / 1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Geoid of Ceres (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');

set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% SHAPE
file = 'CERES_SPC181019_0256.ICQ';
res = 1;

fid = fopen(file, 'r');
fgetl(fid);
data = fscanf(fid, '%f', [3, Inf])';
fclose(fid);

x = data(:, 1);
y = data(:, 2);
z = data(:, 3);

lon = atan2(y, x); 
lon = mod(lon + 180, 360) - 180;
lat = atan2(z, sqrt(x.^2 + y.^2));
lon = rad2deg(lon); 
lat = rad2deg(lat); 
radius = sqrt(x.^2 + y.^2 + z.^2);

latLim = [-89.5 89.5 res];
lonLim = [-179.5 179.5 res];
latGrid = latLim(1):latLim(3):latLim(2);
lonGrid = lonLim(1):lonLim(3):lonLim(2);

[lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
F = scatteredInterpolant(lon, lat, radius, 'natural', 'nearest');
Topography = F(lonMesh, latMesh);

figure;
imagesc(lonGrid, latGrid, Topography);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Shape of Ceres (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');

set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% GRAVITY
file = 'JGDWN_CER18D_SHA.TAB';
res = 1;

fid = fopen(file, 'r');
header = fgetl(fid);
header = textscan(header, '%f%f%f%f%f%f%f%f', 'Delimiter', ',');
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', ',');
fclose(fid);

header = cell2mat(header);
data_matrix = cell2mat(data);
m = data_matrix(:, 1);  % Degree index m
n = data_matrix(:, 2);  % Order index n
Cmn = data_matrix(:, 3);  % Coefficients Cmn
Smn = data_matrix(:, 4);  % Coefficients Smn

SHcoeff = data_matrix(:, 1:4);
c00 = [0., 0., 1., 0.];
SHcoeff = vertcat(c00, SHcoeff);

Radius = header(1) * 1e3;
GM = header(2) * 1e9;
Height = Radius + 10;
SHbounds = [0 header(4)];

latLim = [-89.5 89.5 res];
lonLim = [-179.5 179.5 res];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));

[data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
gravity_acceleration = sqrt(data.vec.R.^2 + data.vec.T.^2 + data.vec.L.^2);

figure;
imagesc(lonGrid, latGrid, gravity_acceleration);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Gravity of Ceres (m/s^2)');
colormap(turbo); 
set(gca, 'YDir', 'normal');

set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% BOUGUER GRAVITY
G = 6.67430e-11;
rho_B = 1800;

heights = Topography .* 1e3 - geoid;
BOUGUER = (2 * pi * G * rho_B) .* heights .*1e5;

figure;
imagesc(lonGrid, latGrid, BOUGUER);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Bouguer Gravity of Ceres (mGal)');
colormap(turbo); 
set(gca, 'YDir', 'normal');

set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% Save the geoid data to a MAT file
save('Data/Geoid_Ceres.mat', 'geoid');

% Save the shape data to a MAT file
save('Data/Shape_Ceres.mat', 'Topography');

% Save the gravity data to a MAT file
save('Data/Gravity_Ceres.mat', 'gravity_acceleration');

% Save the Bouguer gravity data to a MAT file
save('Data/Bouguer_Gravity_Ceres.mat', 'BOUGUER');