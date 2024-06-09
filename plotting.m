clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

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