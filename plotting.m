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
sigma_cr = 0.5;
E_cr = 5e9;
T_e_cr = 2.75e3;
g_avg = 0.27;  

latLim = [-89.5 89.5 1];
lonLim = [-179.5 179.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
[num_rows, num_cols] = size(Lat_matrix);
uniform_matrix = ones(size(Lat_matrix));

% load data
[topo_data] = load("Topo_Ceres.mat");
[SH_data] = load("SHcoeffs_Ceres.mat");
[gravity_data] = load("Gravity_Ceres.mat");
[geoid_data] = load("Geoid_Ceres.mat");

% process data
h_topo_raw = topo_data.Topography;
h_topo = zeros(size(h_topo_raw));
h_topo(:,end/2+1:end) = h_topo_raw(:,1:end/2);
h_topo(:,1:end/2) = h_topo_raw(:,end/2+1:end);

SHcoeff = SH_data.SHcoeff;
SHbounds = [2 18];
[data] = gravityModule(Lat_matrix, Lon_matrix, height, SHbounds, SHcoeff, r_ceres, mu_ceres);
g_obs = data.vec.R;

g_ceres_raw = gravity_data.gravity_acceleration;
g_ceres = zeros(size(g_ceres_raw));
g_ceres(:,end/2+1:end) = g_ceres_raw(:,1:end/2);
g_ceres(:,1:end/2) = g_ceres_raw(:,end/2+1:end);

r_geoid_raw = geoid_data.Geoid;
r_geoid = zeros(size(r_geoid_raw));
r_geoid(:,end/2+1:end) = r_geoid_raw(:,1:end/2);
r_geoid(:,1:end/2) = r_geoid_raw(:,end/2+1:end);

% Initial guesses
rho_c = 1215; % initial crust density in kg/m^3
rho_m = 2429; % initial mantle density in kg/m^3
t_cr = 37.7e3; % initial reference crust thickness in meters

% Airy isostasy
dr = root_airy(rho_c, rho_m, h_topo);

% perform transformation
cs = GSHA(dr, 179);
sc = cs2sc(cs);

n_sh = 1:size(sc,1);
range = 0:250;
phi = flexural_response(range, g_avg, E_cr, T_e_cr, sigma_cr, r_ceres, rho_c, rho_m);

fig = figure('Position', [100, 100, 600, 400], 'Visible', 'off');
plot(range, phi, 'LineWidth', 2);
xlabel('Degree (n)');
ylabel('Flexural Response \Phi(n)');
set(gca, 'YDir', 'normal');
% set(gca, 'XTick', 0:10:150);
set(gca, 'YTick', 0:0.1:1);
% set(gca, 'XTickLabel', 0:10:150);
set(gca, 'YTickLabel', 0:0.1:1);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');
set(gca, 'XScale', 'log'); % Set x-axis to log scale
line([18 18], ylim, 'Color', 'r', 'LineStyle', '--','LineWidth', 2); % Adjust linestyle and color as needed
grid on;

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('Figures/Flexure', '-dpng', '-r300'); 

function phi_n = flexural_response(n, g, E, T_e, sigma, R, rho_c, rho_m)
    
    D = (E * T_e^3) / (12 * (1 - sigma^2));

    phi_n = (1 + (D ./ (g .* (rho_m - rho_c))) .* (2 .* (n + 1) ./ (2*R)).^4).^(-1);

end

function dr = root_airy(rho_c, rho_m, dh)

    dr = (rho_c /  (rho_m - rho_c)) .* dh;

end

% Imports
m1_data = load("Bouguer_Inversion_Thickness_Ceres.mat");
m1 = m1_data.Thickness_matrix_centered;

m2_data_i = load("airy_thicknesses_initial.mat");
m2_i = m2_data_i.t_total_centr;

m2_data_f = load("airy_thicknesses_refined.mat");
m2_f = m2_data_f.t_total_final_centr;

m3_data_i = load("airy_flex_thicknesses_initial.mat");
m3_i = m3_data_i.t_total_centr;

m3_data_f = load("airy_flex_thicknesses_refined.mat");
m3_f = m3_data_f.t_total_final_centr;

lonLim_centered = [-179.5 179.5 1];
lonGrid_centered = lonLim_centered(1):lonLim_centered(3):lonLim_centered(2);
latLim = [-89.5 89.5 1];
latGrid = latLim(1):latLim(3):latLim(2);

caxis_boundary = [22, 48];

fig = figure('Position', [100, 100, 800, 400], 'Visible', 'off');
imagesc(lonGrid_centered, latGrid, m1 ./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Model 1 Crustal Thickness (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');
caxis(caxis_boundary); 

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('Figures/Model1_Crustal_Thickness', '-dpng', '-r300'); 


fig = figure('Position', [100, 100, 800, 400], 'Visible', 'off');
imagesc(lonGrid_centered, latGrid, m2_i ./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Model 2 Crustal Thickness (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');
caxis(caxis_boundary); 

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('Figures/Model2_Crustal_Thickness', '-dpng', '-r300'); 

fig = figure('Position', [100, 100, 800, 400], 'Visible', 'off');
imagesc(lonGrid_centered, latGrid, m3_i ./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Model 3 Crustal Thickness (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');
caxis(caxis_boundary); 

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('Figures/Model3_Crustal_Thickness', '-dpng', '-r300'); 

