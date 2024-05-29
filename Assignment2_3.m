HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);
Model = struct();

Model.number_of_layers = 1;
Model.name = 'Ceres_model';

% Additional variables
Model.GM = 6.2629e+10;
Model.Re = 470000;
Model.geoid = 'none';
Model.nmax = 18;     

latLim = [-89.5 89.5 1];
lonLim = [-179.5 179.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
Uniform_matrix = ones(size(Lat_matrix));

Initial_height = 40;
Initial_density = 1.4;
[top_bound] = matrix2gmt(Uniform_matrix * 0, Lon_matrix, Lat_matrix);
[bottom_bound] = matrix2gmt(Uniform_matrix * Initial_height, Lon_matrix, Lat_matrix);
[top_density] = matrix2gmt(Uniform_matrix * Initial_density, Lon_matrix, Lat_matrix);
[bottom_density] = matrix2gmt(Uniform_matrix * Initial_density, Lon_matrix, Lat_matrix);

Model.l1.bound = [top_bound];
Model.l1.dens  = [top_density];

Model.l2.bound = [bottom_bound];
%Model.l2.dens  = [bottom_density];

save([HOME '/Data/' Model.name '.mat'],'Model')
[SHcoeff] = model_SH_analysis(Model);
Radius = Model.Re;
GM = Model.GM;
Height = Radius + 10;
SHbounds = [0 Model.nmax];

[data] = gravityModule(Lat_matrix,Lon_matrix,Height,SHbounds,SHcoeff,Radius,GM);
gravity_acceleration = sqrt(data.vec.R.^2 + data.vec.T.^2 + data.vec.L.^2);
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
stop


G = 6.67430e-11;

rho_B = 1800;

heights = Topography .* 1e3 - geoid; 
BOUGUER = (2 * pi * G * rho_B) .* heights .* 1e5; 

%initial_crustal_density = 2700; % Initial crustal density in kg/m^3

max_itr = 100;
tol = 1e-5;
scaling_factor = 0.01;

density_deviation = zeros(size(heights));

for iteration = 1:max_itr
    gravity_model = (2 * pi * G * (rho_B + density_deviation)) .* heights .* 1e5;
    
    gravity_residual = BOUGUER - gravity_model;
    
    if max(abs(gravity_residual(:))) < tol
        fprintf('Convergence reached after %d iterations\n', iteration);
        break;
    end
    
    density_deviation = density_deviation + scaling_factor * gravity_residual ./ (2 * pi * G * heights * 1e5);
    
    if mod(iteration, 10) == 0
        figure;
        imagesc(lonGrid, latGrid, density_deviation);
        colorbar;
        xlabel('Longitude (°)');
        ylabel('Latitude (°)');
        title(['Crustal Density Deviations after ', num2str(iteration), ' Iterations']);
        colormap(turbo); 
        set(gca, 'YDir', 'normal');
        set(gca, 'XTick', -180:45:180);
        set(gca, 'YTick', -90:45:90);
        set(gca, 'XTickLabel', -180:45:180);
        set(gca, 'YTickLabel', -90:45:90);
        set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
        set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
        set(gca, 'TickDir', 'out');
    end
end

figure;
imagesc(lonGrid, latGrid, density_deviation);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Final Crustal Density Deviations (kg/m^3)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');
