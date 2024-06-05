clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

G = 6.67430e-11;

Model = struct();

Model.number_of_layers = 2;
Model.name = 'Ceres_model';

% Additional variables
Model.GM = 6.2629e+10;
Model.Re = 470000;
Model.geoid = 'none';
Model.nmax = 18;     

latLim = [-89.5 89.5 1];
lonLim = [0.5 359.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
Uniform_matrix = ones(size(Lat_matrix));

% Crust layer
Crust_thickness = 40;
Crust_density = 1.25;
[crust_bound] = matrix2gmt(Uniform_matrix .* 0, Lon_matrix, Lat_matrix);
[crust_density] = matrix2gmt(Uniform_matrix .* Crust_density, Lon_matrix, Lat_matrix);
Model.l1.bound = [crust_bound];
Model.l1.dens  = [crust_density];

% Mantle layer
Mantle_density = 2.434;
[mantle_bound] = matrix2gmt(Uniform_matrix .* - Crust_thickness, Lon_matrix, Lat_matrix);
[mantle_density] = matrix2gmt(Uniform_matrix .* Mantle_density, Lon_matrix, Lat_matrix);
Model.l2.bound = [mantle_bound];
Model.l2.dens  = [mantle_density];

% Bottom layer
[bottom_bound] = matrix2gmt(Uniform_matrix .* -Model.Re/1e3, Lon_matrix, Lat_matrix);
Model.l3.bound = [bottom_bound];
%Model.l2.dens  = [bottom_density];

save([HOME '/Data/' Model.name '.mat'],'Model')
[SHcoeff] = model_SH_analysis(Model);
Radius = Model.Re;
GM = Model.GM;
Height = Radius + 0;
SHbounds = [0 Model.nmax];

[data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
gravity_acceleration = sqrt(data.vec.X.^2 + data.vec.Y.^2 + data.vec.Z.^2);
figure;
imagesc(lonGrid, latGrid, gravity_acceleration);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Gravity of Ceres (m/s^2)');
colormap(turbo); 
set(gca, 'YDir', 'normal');

set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% % % % % % % % % % % %