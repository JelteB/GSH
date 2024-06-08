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

function X = rands(sz,sig)
        % BASE CASE
        if nargin == 0
            X = rand;
            return
        end
        % PARSE SIZE OF ARRAY
        if any(sz ~= fix(sz))
            error('Size inputs must be integers');
        end
        if numel(sz) == 1
            sz = sz*[1 1];
        end
        % DEFAULT GAUSSIAN WIDTH (pixels)
        if nargin == 1
            sig = ceil(0.1*max(sz));
        end
        % GENERATE RANDOM NUMBERS
        X = rand(sz);
        % CALCULATE GAUSSIAN KERNEL 
        K = gkern(sz,sig);
        % SMOOTH RANDOM NUMBERS WITH FFT
        X = real(ifftn(fftn(X) .* fftn(K,sz)));
        % SCALE X
        X = X - min(X(:));
        X = X / max(X(:));
    end
    function K = gkern(sz,sigma)
    % GKERN     Calculates Gaussian kernel of size "sz" with width "sig"
        % CALCULATE NUMBER OF DIMENSIONS
        N = numel(sz);
        % CALCULATE GRID CENTER
        cidx = ceil(sz / 2);
        
        % CALCULATE ND GRID
        inds = cell(1, N);
        for i = 1 : N
            inds{i} = 1:sz(i);
        end
        [grid{1:N}] = ndgrid(inds{:});
        % CALCULATE RADIAL DISTANCE 
        % SQUARED FROM CENTER
        RSQ = 0;
        for i = 1 : N
            RSQ = RSQ + (grid{i} - cidx(i)).^2;
        end
        % CALCULATE GAUSSIAN KERNEL
        K = exp(- RSQ / (2*sigma^2)); 
        K = K / sum(K(:));
    end

Smooth_noise = rands([ 180 360 ], 2);
Topography = 10e3 .* (Smooth_noise - mean(Smooth_noise));

Smooth_noise2 = rands([ 180 360 ], 2);
Boundary_matrix = synthetic_crust_thickness .* Uniform_matrix + 1e3 .* (Smooth_noise2 - mean(Smooth_noise2));

% Compute synthetic gravity data using forward model
V_Model_synthetic = segment_2layer_model(Topography, Boundary_matrix, -150e3, ...
    synthetic_crust_density, synthetic_mantle_density, 20e3, Model);
[data_synthetic] = model_SH_synthesis(lonLim, latLim, 0., SHbounds, V_Model_synthetic, Model);
synthetic_gravity_matrix = data_synthetic.vec.R;

% Initialize the inversion with a different boundary matrix
initial_guess_thickness = -40e3;
Smooth_noise3 = rands([ 180 360 ], 2);
Boundary_matrix_inverted = initial_guess_thickness .* Uniform_matrix + 4e3 .* (Smooth_noise3 - mean(Smooth_noise3));

% % % ITERATION
max_itr = 150;
tolerance = 1e-5;
best_value = 1;

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

    if abs(mean(residual_matrix, 'all')) < best_value
        best_boundary_matrix = Boundary_matrix_inverted;
        best_value = abs(mean(residual_matrix, 'all'));
        best_iter = iter;
        disp(['Best iteration: ' num2str(iter)])
    end

end

Thickness_matrix = Topography - Boundary_matrix;
Thickness_matrix_inv = Topography - best_boundary_matrix;

% Evaluate results
disp(['True Thickness mean: ', num2str(mean(Thickness_matrix(:)))]);
disp(['Inverted  Thickness Mean: ', num2str(mean(Thickness_matrix_inv(:)))]);
disp(['Best Residual: ', num2str(best_value)]);

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