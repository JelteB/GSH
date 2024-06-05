clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

G = 6.67430e-11;

latLim = [-89.5 89.5 1];
lonLim = [-179.5 179.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
[num_rows, num_cols] = size(Lat_matrix);

[observations] = load("Bouguer_Heights_Ceres.mat");
[observations_bouger] = load("Bouguer_Anomaly_Ceres.mat");

topo_data = observations.heights;
BA_obs = observations_bouger.Bouguer_anom;

max_itr = 25;
tol = 1e-9;

% Initial guesses
rho_c = 1200; % initial crust density in kg/m^3
rho_m = 2400; % initial mantle density in kg/m^3
h_c0 = 40e3; % initial reference crust thickness in meters

% Flatten the data for optimization
BA_obs_flat = BA_obs(:);
topo_data_flat = topo_data(:);

% Perform optimization using least squares
for iter = 1:max_itr
    % Calculate the change in crustal thickness
    delta_h_c = BA_obs_flat / (2 * pi * G * (rho_m - rho_c));
    % Calculate the residuals
    BA_model = 2 * pi * G * (rho_m - rho_c) * delta_h_c + 2 * pi * G * rho_c * topo_data_flat;
    residuals = BA_obs_flat - BA_model;
    
    % Compute the design matrix (H)
    H = 2 * pi * G * [ -delta_h_c + topo_data_flat, delta_h_c, -(rho_m - rho_c) * ones(length(BA_obs_flat), 1)];
    
    % Compute the parameter updates using the normal equations
    params_update = (H' * H) \ (H' * residuals);
    
    % Update the parameters
    rho_c = rho_c + params_update(1);
    rho_m = rho_m + params_update(2);
    h_c0 = h_c0 + params_update(3);
    
    % Compute the cost (sum of squared residuals)
    cost = sum(residuals.^2);
    
    % Check for convergence
    if sqrt(cost) < tol
        break;
    end
    
    % Display iteration details
    fprintf('Iteration %d: Cost = %f, rho_c = %f, rho_m = %f, h_c0 = %f\n', ...
        iter, cost, rho_c, rho_m, h_c0);
end

% Compute final crustal thickness map using optimized parameters
delta_h_c_opt_flat = BA_obs_flat / (2 * pi * G * (rho_m - rho_c));
h_c_opt_flat = h_c0 + delta_h_c_opt_flat;

% Reshape the results back to the original 2D shape
delta_h_c_opt = reshape(delta_h_c_opt_flat, num_rows, num_cols);
h_c_opt = reshape(h_c_opt_flat, num_rows, num_cols);

str_c = ['Crustal density: ', num2str(rho_c), ' [kg/m³]'];
str_m = ['Mantle density: ', num2str(rho_m), ' [kg/m³]'];

disp(str_c);
disp(str_m);

figure;
imagesc(lonGrid, latGrid, delta_h_c_opt./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Crustal Thickness Deviations (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

figure;
imagesc(lonGrid, latGrid, h_c_opt./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Crustal Thickness (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

%% Define functions

