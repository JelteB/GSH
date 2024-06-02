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
[numRows, numCols] = size(Lat_matrix);
Uniform_matrix = ones(size(Lat_matrix));

Average_thickness = 40e3;
Rho_B = 1800;

Average_matrix = Uniform_matrix .* Average_thickness;
Thickness_matrix = Average_matrix .* 1;
[observations] = load("Bouguer_Heights_Ceres.mat");
Heights = observations.heights;

max_itr = 10;
tol = 1e-9;

for itr = 1:max_itr
    Bouguer_obs  = (2 * pi * G * Rho_B) .* Heights;

    % Calculate the deviation from the average thickness
    Deviation_matrix = Thickness_matrix - Average_matrix;
    disp('deviation')
    disp(Deviation_matrix(1,1))

    % Forward calculation of the Bouguer anomaly using current thickness model
    Bouguer_calc = (2 * pi * G * Rho_B) .* Deviation_matrix;
    disp('bouguer calc')
    disp(Bouguer_calc(1,1))

    % Compute residuals between observed and calculated Bouguer anomalies
    Residual_matrix = Bouguer_obs - Bouguer_calc;
    disp('bouguer obs')
    disp(Bouguer_obs(1,1))
    disp('res')
    disp(Residual_matrix(1,1))

    % Correction to the thickness model
    Correction_matrix = Residual_matrix / (2 * pi * G * Rho_B);
    disp('corr')
    disp(Correction_matrix(1,1))

    % Update thickness model
    Thickness_matrix = Thickness_matrix + Correction_matrix;
    disp('thickness')
    disp(Thickness_matrix(1,1))

    % Check for convergence
    if max(abs(Correction_matrix), [], 'all') < tol
        disp(['Converged in ', num2str(itr), ' iterations']);
        break;
    end
end
% for itr = 1:max_itr
%     disp(Bouguer_obs(50,73))
%     Deviation_matrix = Thickness_matrix - Average_matrix;
%     disp(Deviation_matrix(50,73))
%     Residual_matrix = Bouguer_obs - (2 * pi * G * rho_B .* Deviation_matrix);
%     test = (2 * pi * G * rho_B .* Deviation_matrix);
%     disp(test(50,73));
%     disp(Residual_matrix(50,73))
%     Correction_matrix = Residual_matrix / (2 * pi * G * rho_B);
%     disp(Correction_matrix(50,73))
%     Thickness_matrix = Thickness_matrix + Correction_matrix;
%     disp(Thickness_matrix(50,73))

%     if max(abs(Correction_matrix), [], 'all') < tol
%         disp(['Converged in ', num2str(itr), ' iterations']);
%         break;
%     end
%     disp('===========================')
%     stop
% end

Deviation_matrix = Thickness_matrix - Average_matrix;

figure;
imagesc(lonGrid, latGrid, Thickness_matrix./1e3);
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
stop

for itr = 1:max_itr

    Bouguer_itr = (2 * pi * G * rho_B) .* height_matrix;
    
    gravity_residual = Bouguer_obs - Bouguer_itr;

    new_height_matrix = height_matrix + 1 * (gravity_residual / (2 * pi * G * rho_B));

    if max(abs(new_height_matrix - height_matrix), [], 'all') < tol
        disp(['Converged in ', num2str(itr), ' iterations']);
        break;
    end

    height_matrix = new_height_matrix;

end

disp(height_matrix(:,1));

figure;
imagesc(lonGrid, latGrid, height_matrix./1e3);
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

% heights = Topography .* 1e3 - geoid; 
% BOUGUER = (2 * pi * G * rho_B) .* heights .* 1e5; 

% %initial_crustal_density = 2700; % Initial crustal density in kg/m^3

% max_itr = 100;
% tol = 1e-5;
% scaling_factor = 0.01;

% density_deviation = zeros(size(heights));

% for iteration = 1:max_itr
%     gravity_model = (2 * pi * G * (rho_B + density_deviation)) .* heights .* 1e5;
    
%     gravity_residual = BOUGUER - gravity_model;
    
%     if max(abs(gravity_residual(:))) < tol
%         fprintf('Convergence reached after %d iterations\n', iteration);
%         break;
%     end
    
%     density_deviation = density_deviation + scaling_factor * gravity_residual ./ (2 * pi * G * heights * 1e5);
    
%     if mod(iteration, 10) == 0
%         figure;
%         imagesc(lonGrid, latGrid, density_deviation);
%         colorbar;
%         xlabel('Longitude (°)');
%         ylabel('Latitude (°)');
%         title(['Crustal Density Deviations after ', num2str(iteration), ' Iterations']);
%         colormap(turbo); 
%         set(gca, 'YDir', 'normal');
%         set(gca, 'XTick', -180:45:180);
%         set(gca, 'YTick', -90:45:90);
%         set(gca, 'XTickLabel', -180:45:180);
%         set(gca, 'YTickLabel', -90:45:90);
%         set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
%         set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
%         set(gca, 'TickDir', 'out');
%     end
% end

% figure;
% imagesc(lonGrid, latGrid, density_deviation);
% colorbar;
% xlabel('Longitude (°)');
% ylabel('Latitude (°)');
% title('Final Crustal Density Deviations (kg/m^3)');
% colormap(turbo); 
% set(gca, 'YDir', 'normal');
% set(gca, 'XTick', -180:45:180);
% set(gca, 'YTick', -90:45:90);
% set(gca, 'XTickLabel', -180:45:180);
% set(gca, 'YTickLabel', -90:45:90);
% set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
% set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
% set(gca, 'TickDir', 'out');

% Model = struct();

% Model.number_of_layers = 1;
% Model.name = 'Ceres_model';

% % Additional variables
% Model.GM = 6.2629e+10;
% Model.Re = 470000;
% Model.geoid = 'none';
% Model.nmax = 18;     

% latLim = [-89.5 89.5 1];
% lonLim = [-179.5 179.5 1];
% lonGrid = lonLim(1):lonLim(3):lonLim(2);
% latGrid = latLim(1):latLim(3):latLim(2);
% Lon_matrix = repmat(lonGrid,length(latGrid),1);
% Lat_matrix = repmat(latGrid',1,length(lonGrid));
% Uniform_matrix = ones(size(Lat_matrix));

% Initial_height = 40;
% Initial_density = 1.4;
% [top_bound] = matrix2gmt(Uniform_matrix * 0, Lon_matrix, Lat_matrix);
% [bottom_bound] = matrix2gmt(Uniform_matrix * Initial_height, Lon_matrix, Lat_matrix);
% [top_density] = matrix2gmt(Uniform_matrix * Initial_density, Lon_matrix, Lat_matrix);
% [bottom_density] = matrix2gmt(Uniform_matrix * Initial_density, Lon_matrix, Lat_matrix);

% Model.l1.bound = [top_bound];
% Model.l1.dens  = [top_density];

% Model.l2.bound = [bottom_bound];
% %Model.l2.dens  = [bottom_density];

% save([HOME '/Data/' Model.name '.mat'],'Model')
% [SHcoeff] = model_SH_analysis(Model);
% Radius = Model.Re;
% GM = Model.GM;
% Height = Radius + 10;
% SHbounds = [0 Model.nmax];

% [data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
% gravity_acceleration = sqrt(data.vec.R.^2 + data.vec.T.^2 + data.vec.L.^2);
% figure;
% imagesc(lonGrid, latGrid, gravity_acceleration);
% colorbar;
% xlabel('Longitude (°)');
% ylabel('Latitude (°)');
% title('Gravity of Ceres (m/s^2)');
% colormap(turbo); 
% set(gca, 'YDir', 'normal');

% set(gca, 'XTick', -180:45:180);
% set(gca, 'YTick', -90:45:90);
% set(gca, 'XTickLabel', -180:45:180);
% set(gca, 'YTickLabel', -90:45:90);
% set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
% set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
% set(gca, 'TickDir', 'out');

