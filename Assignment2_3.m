clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Setup
latLim = [-89.5 89.5 1];
lonLim = [0.5 359.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
Uniform_matrix = ones(size(Lat_matrix));
G = 6.67430e-11;
gt = 0;


% Imports
[Topo_data] = load("Topo_Ceres.mat");
Topography_raw = Topo_data.Topography;
Topography = zeros(size(Topography_raw));
Topography(:,end/2+1:end) = Topography_raw(:,1:end/2);
Topography(:,1:end/2) = Topography_raw(:,end/2+1:end);

[SH_data] = load("SHcoeffs_Ceres.mat");
SHcoeff = SH_data.SHcoeff;

% % % OBSERVATION
Radius = 470000;
GM = 6.262905361210000e+10;
Height = Radius + 0.;
SHbounds = [2 18];

[data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
obs_gravity_matrix = data.vec.R;

% % % MODEL
Model = struct();
Model.number_of_layers = 2;
Model.name = 'Ceres_model';
Model.GM = GM;
Model.Re = Radius;
Model.geoid = 'none';
Model.nmax = 18;     

% Crust layer
Crust_thickness = - 37.7e3;
Crust_density = 1215;

% Mantle layer
Mantle_density = 2429;
Boundary_matrix = Uniform_matrix .* Crust_thickness;

% Bottom layer
Bottom_thickness = -175e3;

% Extract data
save([HOME '/Data/' Model.name '.mat'],'Model')

% % % ITERATION
max_itr = 350;
tolerance = 1e-5;
best_value = 1;

for iter = 0:max_itr

    [V_Model] = segment_2layer_model(Topography,Boundary_matrix,Bottom_thickness,Crust_density,Mantle_density,20e3,Model);

    V = V_Model;
    [data] = model_SH_synthesis(lonLim,latLim,0.,SHbounds,V,Model);
    model_gravity_matrix = data.vec.R;

    % residual
    residual_matrix = obs_gravity_matrix - model_gravity_matrix; % m/s2

    % new boundary
    Addition_matrix = 0.2 .* (residual_matrix .* 1e5);

    % Update boundary matrix based on residual
    for i = 1:size(residual_matrix, 1)
        for j = 1:size(residual_matrix, 2)
            % Only update if residual exceeds tolerance
            if abs(residual_matrix(i, j)) > tolerance
                Boundary_matrix(i, j) = Boundary_matrix(i, j) + Addition_matrix(i, j);
            end
        end
    end

    if abs(mean(residual_matrix, 'all')) < best_value
        best_boundary_matrix = Boundary_matrix;
        best_value = abs(mean(residual_matrix, 'all'));
    end

end

Thickness_matrix = Topography - best_boundary_matrix;
[Thickness_matrix_centered] = Europe_centered(Thickness_matrix);

lonLim_centered = [-179.5 179.5 1];
lonGrid_centered = lonLim_centered(1):lonLim_centered(3):lonLim_centered(2);

figure('Position',[100 100 800 400]);
imagesc(lonGrid_centered, latGrid(1+gt:end-gt), Thickness_matrix_centered(1+gt:end-gt,:) ./1e3);
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
caxis_boundary = [2.507429118515038e+01, 4.735139342625794e+01];
caxis(caxis_boundary); % Set color scale range for boundary

% Save Tickness
save('Data/Bouguer_Inversion_Thickness_Ceres.mat', 'Thickness_matrix_centered');