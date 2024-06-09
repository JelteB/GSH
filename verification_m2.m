clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% define parameters
G = 6.67430e-11;
mu_ceres = 6.262905361210000e+10;
r_ceres = 470e3;
height = r_ceres;
gt = 0;
scale_param = 20;

latLim = [-89.5 89.5 1];
lonLim = [-179.5 179.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
[num_rows, num_cols] = size(Lat_matrix);
uniform_matrix = ones(size(Lat_matrix));

% data
rho_c_p = 1200; 
rho_m_p = 2000; 
rho_c_n = 800; 
rho_m_n = 1600; 
t_cr = 20e3; 

% Set up synthetic topographic profiles
amplitude = 0.5e3;
frequency = 0.1; 

[x, y] = meshgrid(1:num_cols, 1:num_rows);

% Create the surface using a combination of sine and cosine functions
dh = 0.5 * amplitude * (sin(frequency * x) + cos(frequency * y));
h_topo_pos = dh + 2 * amplitude;
h_topo_neg = dh - 2 * amplitude; 
h_topo = dh; 

dr_pos = root_airy(rho_c_p, rho_m_p, h_topo_pos);
dr_neg = root_airy(rho_c_n, rho_m_n, h_topo_neg);

t_total_neg = t_cr + h_topo_neg + dr_neg;
t_total_pos = t_cr + h_topo_pos + dr_pos;

%% plots

figure('Units', 'inches', 'Position', [0, 0, 24, 12]);

subplot(2, 2, 1);
imagesc(lonGrid, latGrid, t_total_neg ./ 1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Crustal thickness for negative elevations (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

subplot(2, 2, 2);
imagesc(lonGrid, latGrid, t_total_pos ./ 1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Crustal thickness for positive elevations (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

subplot(2, 2, 3);
imagesc(lonGrid, latGrid, h_topo_neg ./ 1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Negative topographic elevation (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

subplot(2, 2, 4);
imagesc(lonGrid, latGrid, h_topo_pos ./ 1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Positive topographic elevation (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

h = gcf;
h.Position(3) = 24; % Set figure width
h.Position(4) = 12; % Set figure height
h.PaperSize = [24,12]; % Set paper size for printing
h.PaperPositionMode = 'manual';
h.PaperPosition = [0, 0, 24, 12]; % Set paper position


%% set up functions

function dr = root_airy(rho_c, rho_m, dh)

    dr = (rho_c /  (rho_m - rho_c)) .* dh;

end