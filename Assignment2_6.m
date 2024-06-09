clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Imports
obs_data = load("SHcoeffs_Ceres.mat");
Obs = obs_data.SHcoeff;
M1_data = load("M1_SH.mat");
M1 = M1_data.best_SH_coeff;

M2_initial_data = load("M2_SH_initial.mat");
M2_i = M2_initial_data.V_initial;
M2_final_data = load("M2_SH_final.mat");
M2_f = M2_final_data.V_final;

M3_initial_data = load("M3_SH_initial.mat");
M3_i = M3_initial_data.V_initial;
M3_final_data = load("M3_SH_final.mat");
M3_f = M3_final_data.V_final;

% Create a figure with defined size
figure('Position',[100 100 600 400], 'Visible', 'off');

% Hold on to allow multiple plots on the same figure
hold on;

% Plot all models
plot_DV(Obs, 'Observed', 'red', '-');
plot_DV(M1, 'Model 1', 'yellow', '-');
plot_DV(M2_i, 'Model 2', '#7E2F8E', '-');
% plot_DV(M2_f, 'Model 2 (f)', '#7E2F8E', ':');
plot_DV(M3_i, 'Model 3', '#4DBEEE', '-');
% plot_DV(M3_f, 'Model 3 (f)', '#4DBEEE', ':');

% Labels and title
xlabel('Degree (-)');
ylabel('Degree Variance (-)');
title('Degree Variance of Models');
set(gca, 'YScale', 'log');
legend;
grid on;

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('Figures/DV_plot', '-dpng', '-r300'); 

% Function to plot degree variance
function plot_DV(Sh_coeff_matrix, name, color, style)
    [n, DV] = degreeVariance(Sh_coeff_matrix);
    % Plot using line instead of scatter for better visibility
    plot(n, DV, 'DisplayName', name, 'LineWidth', 2, 'LineStyle', style, 'Color', color);
end

