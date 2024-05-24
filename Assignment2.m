clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% TOPOGRAPHY
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

latGrid = linspace(-90, 90, 180 * res); 
lonGrid = linspace(-180, 180, 360 * res + 1); % Adjusted for -180 to 180 range
[lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
F = scatteredInterpolant(lon, lat, radius, 'natural', 'nearest');
radiusGrid = F(lonMesh, latMesh);

figure;
imagesc(lonGrid, latGrid, radiusGrid);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Topography of Ceres (km)');
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
lonLim = [-180 180 res];
lon = lonLim(1):lonLim(3):lonLim(2);
lat = latLim(1):latLim(3):latLim(2);

Lon_matrix = repmat(lon,length(lat),1);
Lat_matrix = repmat(lat',1,length(lon));

[data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
gravity_acceleration = sqrt(data.vec.R.^2 + data.vec.T.^2 + data.vec.L.^2);

figure;
imagesc(lon, lat, gravity_acceleration);
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

disp(size(gravity_acceleration));
disp(size(radiusGrid));

% lmax = max(data_matrix(:, 1));
% field = zeros(lmax + 1, 2*lmax+1);
% mid_row = ceil((2*lmax+1)/ 2);

% for i = 1:size(data_matrix, 1)
%     m_val = m(i);
%     n_val = n(i);
%     Cmn_val = Cmn(i);
%     Smn_val = Smn(i);

%     if n_val == 0
%         field(m_val + 1, mid_row) = Cmn_val;
%     else
%         field(m_val + 1, n_val + mid_row) = Cmn_val;
%         field(m_val + 1, - n_val + mid_row) = Smn_val;
%     end

% end

% lam = [res/2:res:360-res/2]; % lam   [n x 1]   longitude [deg]
% th = [res/2:res:180-res/2]; % th    [m x 1]   co-latitude [deg]      【90-lat】
% test = GSHS(field,lam,th);