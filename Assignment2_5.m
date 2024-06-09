clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% define parameters
G = 6.67430e-11;
mu_ceres = 6.262905361210000e+10;
r_ceres = 470e3;
height = r_ceres;
gt = 0;
scale_param = 20;
sigma_cr = 0.5;
E_cr = 5e9;
T_e_cr = 2.75e3;
g_avg = 0.27;  

latLim = [-89.5 89.5 1];
lonLim = [-179.5 179.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
[num_rows, num_cols] = size(Lat_matrix);
uniform_matrix = ones(size(Lat_matrix));

% load data
[topo_data] = load("Topo_Ceres.mat");
[SH_data] = load("SHcoeffs_Ceres.mat");
[gravity_data] = load("Gravity_Ceres.mat");
[geoid_data] = load("Geoid_Ceres.mat");

% process data
h_topo_raw = topo_data.Topography;
h_topo = zeros(size(h_topo_raw));
h_topo(:,end/2+1:end) = h_topo_raw(:,1:end/2);
h_topo(:,1:end/2) = h_topo_raw(:,end/2+1:end);

SHcoeff = SH_data.SHcoeff;
SHbounds = [2 18];
[data] = gravityModule(Lat_matrix, Lon_matrix, height, SHbounds, SHcoeff, r_ceres, mu_ceres);
g_obs = data.vec.R;

g_ceres_raw = gravity_data.gravity_acceleration;
g_ceres = zeros(size(g_ceres_raw));
g_ceres(:,end/2+1:end) = g_ceres_raw(:,1:end/2);
g_ceres(:,1:end/2) = g_ceres_raw(:,end/2+1:end);

r_geoid_raw = geoid_data.Geoid;
r_geoid = zeros(size(r_geoid_raw));
r_geoid(:,end/2+1:end) = r_geoid_raw(:,1:end/2);
r_geoid(:,1:end/2) = r_geoid_raw(:,end/2+1:end);

% Initial guesses
rho_c = 1215; % initial crust density in kg/m^3
rho_m = 2429; % initial mantle density in kg/m^3
t_cr = 37.7e3; % initial reference crust thickness in meters

% Airy isostasy
dr = root_airy(rho_c, rho_m, h_topo);

% perform transformation
cs = GSHA(dr, 179);
sc = cs2sc(cs);

n_sh = 1:size(sc,1);

phi = flexural_response(n_sh, g_avg, E_cr, T_e_cr, sigma_cr, r_ceres, rho_c, rho_m);

sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*phi';
end

dr_flex = GSHS(sc_flex,lonGrid,90-latGrid,179);


t_total = t_cr + h_topo + dr_flex;

[t_total_centr] = Europe_centered(t_total);

str_c = ['Crustal density: ', num2str(rho_c), ' [kg/m³]'];
str_m = ['Mantle density: ', num2str(rho_m), ' [kg/m³]'];

disp(str_c);
disp(str_m);

% iteration parameters
max_itr = 100;
tol = 1e-5;
t_boundary = - uniform_matrix .* t_cr - dr;
res_mean_prev = uniform_matrix;

for iter = 0:max_itr   

    g_model = gravity_model(rho_c, rho_m, h_topo, t_boundary, r_ceres, mu_ceres);

    % residual
    residual_matrix = g_obs - g_model;
    
    res_max = max(residual_matrix, [], "all");    
    str_mean = ["Max residual: ", num2str(res_max)];
    str_iter = ["Iteration number: ", num2str(iter)];

    disp(str_iter);
    disp(str_mean);

    % scaling
    addition_matrix = 0.2 .* (residual_matrix .* 1e5);

    % Update boundary matrix based on residual
    for i = 1:size(residual_matrix, 1)
        for j = 1:size(residual_matrix, 2)
            % Only update if residual exceeds tolerance
            if abs(residual_matrix(i, j)) > tol
                t_boundary(i, j) = t_boundary(i, j) + addition_matrix(i, j);
            end
        end
    end
end

t_total_final = - t_boundary + h_topo;
t_total_final_centr = Europe_centered(t_total_final);

lonLim_centered = [-179.5 179.5 1];
lonGrid_centered = lonLim_centered(1):lonLim_centered(3):lonLim_centered(2);


figure('Position',[100 100 800 400]);
imagesc(lonGrid_centered, latGrid(1+gt:end-gt), t_total_centr(1+gt:end-gt,:) ./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Initial Crustal Thickness (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

figure('Position',[100 100 800 400]);
imagesc(lonGrid_centered, latGrid(1+gt:end-gt), t_total_final_centr(1+gt:end-gt,:) ./1e3);
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

%% save thicknesses
save('Data/airy_flex_thicknesses_refined.mat', 't_total_final_centr');
save('Data/airy_flex_thicknesses_initial.mat', 't_total_centr');


%% Define functions

function dr = root_airy(rho_c, rho_m, dh)

    dr = (rho_c /  (rho_m - rho_c)) .* dh;

end

function g_model = gravity_model(rho_c, rho_m, h_top, t_b, R, mu)

    % model definition
    Model = struct();    
    Model.number_of_layers = 2;
    Model.name = 'Ceres_model';
    Model.GM = mu;
    Model.Re = R;
    Model.geoid = 'none';
    Model.nmax = 18;     
    
    latLim = [-89.5 89.5 1];
    lonLim = [0.5 359.5 1];

    % layer boundaries
    boundary_matrix = (t_b);
    t_bottom = -175e3;

    % [V] = model_SH_analysis(Model);
    V = segment_2layer_model(h_top, boundary_matrix, t_bottom, rho_c, rho_m, 20e3, Model);
    SHbounds = [2 18];
    height = R;
    [data] = model_SH_synthesis(lonLim,latLim, 0.,SHbounds,V,Model);
    g_model = data.vec.R;

end


function phi_n = flexural_response(n, g, E, T_e, sigma, R, rho_c, rho_m)
    
    D = (E * T_e^3) / (12 * (1 - sigma^2));

    phi_n = (1 + (D ./ (g .* (rho_m - rho_c))) .* (2 .* (n + 1) ./ (2*R)).^4).^(-1);

end