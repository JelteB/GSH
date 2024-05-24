clear;
close all;
clc;

HOME = pwd;

addpath([HOME '/Data']);
addpath([HOME '/Tools']);


% % TOPOGRAPHY
% file = 'CERES_SPC181019_0256.ICQ';
% res = 1;

% fid = fopen(file, 'r');
% fgetl(fid);
% data = fscanf(fid, '%f', [3, Inf])';
% fclose(fid);

% x = data(:, 1);
% y = data(:, 2);
% z = data(:, 3);

% lon = atan2(y, x); 
% lon = mod(lon, 2*pi); 
% lat = atan2(z, sqrt(x.^2 + y.^2));
% lon = rad2deg(lon); 
% lat = rad2deg(lat); 
% radius = sqrt(x.^2 + y.^2 + z.^2);

% latGrid = linspace(-90, 90, 180 * res + 1); 
% lonGrid = linspace(0, 360, 360 * res + 1); 
% [lonMesh, latMesh] = meshgrid(lonGrid, latGrid);
% F = scatteredInterpolant(lon, lat, radius, 'natural', 'nearest');
% radiusGrid = F(lonMesh, latMesh);

% figure;
% imagesc(lonGrid, latGrid, radiusGrid);
% colorbar;
% xlabel('Longitude (°)');
% ylabel('Latitude (°)');
% title('2D Topography of Ceres (Radius)');
% colormap(turbo); 
% set(gca, 'YDir', 'normal');
% axis equal;

% set(gca, 'XTick', 0:45:360);
% set(gca, 'YTick', -90:45:90);
% set(gca, 'XTickLabel', 0:45:360);
% set(gca, 'YTickLabel', -90:45:90);
% set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
% set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
% set(gca, 'TickDir', 'out');

% GRAVITY
file = 'JGDWN_CER18D_SHA.TAB';
res = 1;

fid = fopen(file, 'r');
header = fgetl(fid);
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', ',');
fclose(fid);

data_matrix = cell2mat(data);

m = data_matrix(:, 1);  % Degree index m
n = data_matrix(:, 2);  % Order index n
Cmn = data_matrix(:, 3);  % Coefficients Cmn
Smn = data_matrix(:, 4);  % Coefficients Smn
uncertainty_Cmn = data_matrix(:, 5);  % Uncertainties in Cmn
uncertainty_Smn = data_matrix(:, 6);  % Uncertainties in Smn

table_data = table(m, n, Cmn, Smn, uncertainty_Cmn, uncertainty_Smn);

lmax = max(data_matrix(:, 1));
field = zeros(lmax + 1, 2*lmax+1);
mid_row = ceil((2*lmax+1)/ 2);

for i = 1:size(data_matrix, 1)
    m_val = m(i);
    n_val = n(i);
    Cmn_val = Cmn(i);
    Smn_val = Smn(i);

    if n_val == 0
        field(m_val + 1, mid_row) = Cmn_val;
    else
        field(m_val + 1, n_val + mid_row) = Cmn_val;
        field(m_val + 1, - n_val + mid_row) = Smn_val;
    end

end

lam = [res/2:res:360-res/2]; % lam   [n x 1]   longitude [deg]
th = [res/2:res:180-res/2]; % th    [m x 1]   co-latitude [deg]      【90-lat】
test = GSHS(field,lam,th);

% % Plot the gravity data
% figure;
% imagesc(lam, 90 - th, test);
% colorbar; 
% xlabel('Longitude (°)');
% ylabel('Latitude (°)');
% title('Gravity Map');
% set(gca, 'YDir', 'normal');

% Lat = (linspace(-90, 90, 180 * res + 1));
% Lon = (linspace(0, 360, 360 * res + 1));
V = data_matrix(:, 1:4);
c00 = [0., 0., 1., 0.];
V = vertcat(c00, V);
Re = 470000;
GM = 626290536121;
r = Re + 10;
SHbounds = [0 18];

latLim = [-89.5 89.5 1];
lonLim = [-180 180 1];
lon = lonLim(1):lonLim(3):lonLim(2);
lat = latLim(1):latLim(3):latLim(2);
disp(lat)

Lon = repmat(lon,length(lat),1);
Lat = repmat(lat',1,length(lon));

[data] = gravityModule(Lat,Lon,r,SHbounds,V,Re,GM);
% input: Lat: latitude in degree [matrix]
%        Lon: longitude in degree [marix]
%        r: radial distance of computation surface (r=Re+h) [scalar]
%        SHbounds: [minSH max SH] example [0 5] or [2 7]
%        V : SH coefficients [same size as setleg files]
%        V format:  degree; order; Cnm; Snm
%        Re: radius of the Earth [meters]
%        GM: gravitational parameter of Earth 
%
% output: data: data structure of the gravity field
% V:                Spherical harmonic coefficients (4pi-normalization)
%     format:           [l,m,Clm,Slm]     l = 0,1,2,3,... m = 0,0,0,0,... etc.

% latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
% lonLim =    [-180 180 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
% height =    10.0; % height of computation above spheroid
% SHbounds =  [0 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

% [data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
gravity_acceleration = sqrt(data.vec.R.^2 + data.vec.T.^2 + data.vec.L.^2);

figure;
imagesc(lon, lat, data.vec.R);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Gravity Acceleration Map (m/s^2)');