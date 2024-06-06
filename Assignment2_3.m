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

% [Geoid_data] = load("Geoid_Ceres.mat");
% Geoid = Geoid_data.Geoid;
% [Bouguer_data] = load("Bouguer_Anomaly_Ceres.mat");
% Bouguer = Bouguer_data.Bouguer_anom;

[SH_data] = load("SHcoeffs_Ceres.mat");
SHcoeff = SH_data.SHcoeff;

% % % OBSERVATION
Radius = 470000;
GM = 6.262905361210000e+10;
Height = Radius + 0.;
SHbounds = [2 18];

[data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
obs_gravity_matrix = data.vec.R;%sqrt(data.vec.R .^2 + data.vec.T .^2 + data.vec.L .^2);

figure;
imagesc(lonGrid, latGrid(1+gt:end-gt), obs_gravity_matrix(1+gt:end-gt,:));
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('OBS Gravity (m/s^2)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', 0:45:360);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', 0:45:360);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% % % MODEL
Model = struct();
Model.number_of_layers = 2;
Model.name = 'Ceres_model';
Model.GM = GM;
Model.Re = Radius;
Model.geoid = 'none';
Model.nmax = 18;     

% Crust layer
Crust_thickness = - 40e3;
Crust_density = 1250;
% [crust_bound] = matrix2gmt(Topography ./1e3, Lon_matrix, Lat_matrix);
% [crust_density] = matrix2gmt(Uniform_matrix .* Crust_density, Lon_matrix, Lat_matrix);
% Model.l1.bound = [crust_bound];
% Model.l1.dens  = [crust_density];

% Mantle layer
Mantle_density = 2434;
Boundary_matrix = Uniform_matrix .* Crust_thickness;
% [mantle_bound] = matrix2gmt(Boundary_matrix, Lon_matrix, Lat_matrix);
% [mantle_density] = matrix2gmt(Uniform_matrix .* Mantle_density, Lon_matrix, Lat_matrix);
% Model.l2.bound = [mantle_bound];
% Model.l2.dens  = [mantle_density];

% Bottom layer
Bottom_thickness = -175e3;
% [bottom_bound] = matrix2gmt(Uniform_matrix .* Bottom_thickness, Lon_matrix, Lat_matrix);
% Model.l3.bound = [bottom_bound];

% Extract data
%save([HOME '/Data/' Model.name '.mat'],'Model')

% % % ITERATION
max_itr = 75;
tolerance = 1e-5;

for iter = 0:max_itr

    [V_Model] = segment_2layer_model(Topography,Boundary_matrix,Bottom_thickness,Crust_density,Mantle_density,20e3,Model);
    % input:            - top_bound             nxm matrix[in meters]
    %                   - middle_bound          nxm matrix[in meters]
    %                   - bottom_bound          scalar[in meters]
    %                   - dens1                 nxm matrix, or scalar[in kg/m3]
    %                   - dens2                 nxm matrix, or scalar[in kg/m3]
    %                   - thickness_segment     scalar[in meters]
    %                   - Model                 structure for Model.GM and Model.Re
    %
    % output:           - V_Model               nX4 vector containing SH coefficients

    V = V_Model;
    [data] = model_SH_synthesis(lonLim,latLim,0.,SHbounds,V,Model);
    model_gravity_matrix = data.vec.R;%sqrt(data.vec.R .^2 + data.vec.T .^2 + data.vec.L .^2);

    % residual
    residual_matrix = obs_gravity_matrix - model_gravity_matrix; % m/s2

    % % Check
    % if max(max(residual_matrix)) < 1e-6
    %     disp('broken')
    %     break;
    % end

    % new boundary
    Addition_matrix = 0.2 .* (residual_matrix .* 1e5);
    %Boundary_matrix = Boundary_matrix + Addition_matrix;

    % Update boundary matrix based on residual
    for i = 1:size(residual_matrix, 1)
        for j = 1:size(residual_matrix, 2)
            % Only update if residual exceeds tolerance
            if abs(residual_matrix(i, j)) > tolerance
                Boundary_matrix(i, j) = Boundary_matrix(i, j) + Addition_matrix(i, j);
            end
        end
    end

    % plot
    if mod(iter, 25) == 0

        %disp(['iteration ' num2str(iter)])
        disp(['gravity obs ' num2str(obs_gravity_matrix(94,72)) ' ... ' num2str(obs_gravity_matrix(134,71))])
        disp(['gravity model ' num2str(model_gravity_matrix(94,72)) ' ... ' num2str(model_gravity_matrix(134,71))])
        disp(['residual ' num2str(residual_matrix(94,72)) ' ... ' num2str(residual_matrix(134,71))])
        disp(['boundary O ' (num2str(Boundary_matrix(94,72) - Addition_matrix(94,72))) ' ... ' (num2str(Boundary_matrix(134,71) - Addition_matrix(134,71)))])
        disp(['addition ' num2str(Addition_matrix(94,72)) ' ... ' num2str(Addition_matrix(134,71))])
        disp(['boundary N ' num2str(Boundary_matrix(94,72)) ' ... ' num2str(Boundary_matrix(134,71))])
        disp(['res mean ' num2str(mean(residual_matrix, 'all'))])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        figure('Position', [100 100 1200 800]);

        % Subplot 1: Model Gravity
        subplot(2, 2, 1);
        imagesc(lonGrid, latGrid(1+gt:end-gt), model_gravity_matrix(1+gt:end-gt,:));
        colorbar;
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title(['Model Gravity (m/s^2) - Iteration ', num2str(iter)]);
        colormap(turbo); 
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', 0:45:360);
        set(gca, 'YTick', -90:45:90);
        set(gca, 'XTickLabel', 0:45:360);
        set(gca, 'YTickLabel', -90:45:90);
        set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
        set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
        set(gca, 'TickDir', 'out');

        % Subplot 2: Observed Gravity
        subplot(2, 2, 2);
        imagesc(lonGrid, latGrid(1+gt:end-gt), Addition_matrix(1+gt:end-gt,:));
        colorbar;
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title(['Addition (m) - Iteration ', num2str(iter)]);
        colormap(turbo); 
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', 0:45:360);
        set(gca, 'YTick', -90:45:90);
        set(gca, 'XTickLabel', 0:45:360);
        set(gca, 'YTickLabel', -90:45:90);
        set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
        set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
        set(gca, 'TickDir', 'out');

        % Subplot 3: Residual Matrix
        subplot(2, 2, 3);
        imagesc(lonGrid, latGrid(1+gt:end-gt), residual_matrix(1+gt:end-gt,:));
        colorbar;
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title(['Residual Matrix (m/s^2) - Iteration ', num2str(iter)]);
        colormap(turbo); 
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', 0:45:360);
        set(gca, 'YTick', -90:45:90);
        set(gca, 'XTickLabel', 0:45:360);
        set(gca, 'YTickLabel', -90:45:90);
        set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
        set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
        set(gca, 'TickDir', 'out');

        % Subplot 4: Boundary Matrix
        subplot(2, 2, 4);
        imagesc(lonGrid, latGrid(1+gt:end-gt), Boundary_matrix(1+gt:end-gt,:));
        colorbar;
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title(['Boundary Matrix (m) - Iteration ', num2str(iter)]);
        colormap(turbo); 
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', 0:45:360);
        set(gca, 'YTick', -90:45:90);
        set(gca, 'XTickLabel', 0:45:360);
        set(gca, 'YTickLabel', -90:45:90);
        set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
        set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
        set(gca, 'TickDir', 'out');
    end

end
    
% figure('Position',[100 100 800 400]);
% imagesc(lonGrid, latGrid, model_gravity_matrix);
% colorbar;
% xlabel('Longitude (°)');
% ylabel('Latitude (°)');
% title('Model Gravity of Ceres (m/s^2)');
% colormap(turbo); 
% set(gca, 'YDir', 'normal');

% set(gca, 'XTick', 0:45:360);
% set(gca, 'YTick', -90:45:90);
% set(gca, 'XTickLabel', 0:45:360);
% set(gca, 'YTickLabel', -90:45:90);
% set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
% set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
% set(gca, 'TickDir', 'out');

% figure('Position',[100 100 800 400]);
% imagesc(lonGrid, latGrid, gravity_acceleration);
% colorbar;
% xlabel('Longitude (°)');
% ylabel('Latitude (°)');
% title('Observed Gravity of Ceres (m/s^2)');
% colormap(turbo); 
% set(gca, 'YDir', 'normal');

% set(gca, 'XTick', 0:45:360);
% set(gca, 'YTick', -90:45:90);
% set(gca, 'XTickLabel', 0:45:360);
% set(gca, 'YTickLabel', -90:45:90);
% set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
% set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
% set(gca, 'TickDir', 'out');

Thickness_matrix = Topography - Boundary_matrix;
[Thickness_matrix_centered] = Europe_centered(Thickness_matrix);

lonLim_centered = [-179.5 179.5 res];
lonGrid_centered = lonLim(1):lonLim(3):lonLim(2);

figure('Position',[100 100 800 400]);
imagesc(lonGrid_centered, latGrid(1+gt:end-gt), Thickness_matrix_centered(1+gt:end-gt,:));
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Thickness (m)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');