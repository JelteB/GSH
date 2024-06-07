clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% % % OBSERVATION
Radius = 470000;
GM = 6.262905361210000e+10;
Height = Radius + 0.;
SHbounds = [2 18];

% % % MODEL
Model = struct();
Model.number_of_layers = 2;
Model.name = 'Ceres_model';
Model.GM = GM;
Model.Re = Radius;
Model.geoid = 'none';
Model.nmax = 18;   

% Define synthetic model parameters
synthetic_crust_thickness = -40e3; % 40 km
synthetic_crust_density = 2700; % kg/m^3
synthetic_mantle_density = 3300; % kg/m^3

% Generate synthetic topography and boundary
latLim = [-89.5 89.5 1];
lonLim = [0.5 359.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
Uniform_matrix = ones(size(Lat_matrix));

Gaussian_noise = randn(size(Uniform_matrix));
min_height = -5000;
max_height = 10000;
Gaussian_noise = (Gaussian_noise - min(Gaussian_noise(:))) / (max(Gaussian_noise(:)) - min(Gaussian_noise(:))) ...
    * (max_height - min_height) + min_height;
Topography = Gaussian_noise;

Boundary_matrix = synthetic_crust_thickness .* Uniform_matrix + 2 .* - Topography;

% Compute synthetic gravity data using forward model
V_Model_synthetic = segment_2layer_model(Topography, Boundary_matrix, -150e3, ...
    synthetic_crust_density, synthetic_mantle_density, 20e3, Model);
[data_synthetic] = model_SH_synthesis(lonLim, latLim, 0., SHbounds, V_Model_synthetic, Model);
synthetic_gravity_matrix = data_synthetic.vec.R;

% Initialize the inversion with a different boundary matrix
initial_guess_thickness = -40e3;
Gaussian_noise = randn(size(Uniform_matrix));
min_height = -10000;
max_height = 10000; 
Gaussian_noise = (Gaussian_noise - min(Gaussian_noise(:))) / (max(Gaussian_noise(:)) - min(Gaussian_noise(:))) ...
    * (max_height - min_height) + min_height;
Boundary_matrix_inverted = initial_guess_thickness .* Uniform_matrix + Gaussian_noise;

% % % ITERATION
max_itr = 50;
tolerance = 1e-5;
for iter = 0:max_itr
    % Update model
    [V_Model] = segment_2layer_model(Topography, Boundary_matrix_inverted, -150e3, ...
        synthetic_crust_density, synthetic_mantle_density, 20e3, Model);
    [data] = model_SH_synthesis(lonLim, latLim, 0., SHbounds, V_Model, Model);
    model_gravity_matrix = data.vec.R;

    % Compute residual
    residual_matrix = synthetic_gravity_matrix - model_gravity_matrix;

    % Print initial residual norm
    if iter == 0
        disp(['Initial Residual Norm: ', num2str(norm(residual_matrix))]);
    end

    % Update boundary matrix
    Addition_matrix = 0.2 .* (residual_matrix .* 1e5);
    for i = 1:size(residual_matrix, 1)
        for j = 1:size(residual_matrix, 2)
            if abs(residual_matrix(i, j)) > tolerance
                Boundary_matrix_inverted(i, j) = Boundary_matrix_inverted(i, j) + Addition_matrix(i, j);
            end
        end
    end

    % Check for convergence
    max_residual = max(max(abs(residual_matrix)));
    disp(['Iteration: ', num2str(iter), ', Max Residual: ', num2str(max_residual)]);
    
    if max_residual < tolerance
        disp(['Converged at iteration: ', num2str(iter)]);
        break;
    end
end



Thickness_matrix = Topography - Boundary_matrix;
Thickness_matrix_inv = Topography - Boundary_matrix_inverted;

% Evaluate results
disp(['True Thickness mean: ', num2str(mean(Thickness_matrix(:)))]);
disp(['Inverted  Thickness Mean: ', num2str(mean(Thickness_matrix_inv(:)))]);
disp(['Final Residual Norm: ', num2str(norm(residual_matrix))]);

% Plotting results
max_boundary = max(max(Thickness_matrix./1e3));
min_boundary = min(min(Thickness_matrix./1e3));
max_gravity = max(max(synthetic_gravity_matrix.*1e5));
min_gravity = min(min(synthetic_gravity_matrix.*1e5));

caxis_boundary = [min_boundary, max_boundary];
caxis_gravity = [min_gravity, max_gravity];

figure('Units', 'inches', 'Position', [0, 0, 24, 12]);

% Plot synthetic crust thickness
subplot(2, 2, 1);
imagesc(lonGrid, latGrid, Thickness_matrix./1e3);
colorbar;
caxis(caxis_boundary); % Set color scale range for boundary
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Synthetic Crust Thickness (km)');
set(gca, 'YDir', 'normal');

% Plot synthetic gravity anomaly
subplot(2, 2, 2);
imagesc(lonGrid, latGrid, synthetic_gravity_matrix.*1e5);
colorbar;
caxis(caxis_gravity); % Set color scale range for gravity
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Synthetic Gravity Anomaly (mGal)');
set(gca, 'YDir', 'normal');

% Plot derived crust thickness
subplot(2, 2, 3);
imagesc(lonGrid, latGrid, Thickness_matrix_inv./1e3);
colorbar;
caxis(caxis_boundary); % Set color scale range for boundary
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Derived Crust Thickness (km)');
set(gca, 'YDir', 'normal');

% Plot derived gravity anomaly
subplot(2, 2, 4);
imagesc(lonGrid, latGrid, model_gravity_matrix.*1e5);
colorbar;
caxis(caxis_gravity); % Set color scale range for gravity
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Derived Gravity Anomaly (mGal)');
set(gca, 'YDir', 'normal');

% Adjust subplot spacing
h = gcf;
h.Position(3) = 24; % Set figure width
h.Position(4) = 12; % Set figure height
h.PaperSize = [24,12]; % Set paper size for printing
h.PaperPositionMode = 'manual';
h.PaperPosition = [0, 0, 24, 12]; % Set paper position