clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Initial guess for Te
Te0 = 2.75e3;

% Optimization options
% opt= options('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'TolX', 1e-12, 'TolFun', 1e-12);

% Lower and upper bounds for Te
lb = 0.5e3;
ub = 5e3;

% Run the optimization
objective = @(Te) residual_func(Te);
[Te_opt, fval] = fminbnd(objective, lb, ub);

% Display the optimal Te
fprintf('Optimal Te: %.2f m\n', Te_opt);

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


function res = residual_func(T_e)
    % define parameters
    G = 6.67430e-11;
    mu_ceres = 6.262905361210000e+10;
    r_ceres = 470e3;
    height = r_ceres;
    gt = 0;
    scale_param = 20;
    sigma_cr = 0.5;
    E_cr = 5e9;
    g_avg = 0.27;  

    Model = struct();    
    Model.number_of_layers = 2;
    Model.name = 'Ceres_model';
    Model.GM = mu_ceres;
    Model.Re = r_ceres;
    Model.geoid = 'none';
    Model.nmax = 18; 

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
    [obs_data] = load("SHcoeffs_Ceres.mat");
    
    % process data
    h_topo_raw = topo_data.Topography;
    h_topo = zeros(size(h_topo_raw));
    h_topo(:,end/2+1:end) = h_topo_raw(:,1:end/2);
    h_topo(:,1:end/2) = h_topo_raw(:,end/2+1:end);
    V_obs = obs_data.SHcoeff;

    SHcoeff = SH_data.SHcoeff;
    SHbounds = [2 18];
    [data] = gravityModule(Lat_matrix, Lon_matrix, height, SHbounds, SHcoeff, r_ceres, mu_ceres);
    g_obs = data.vec.R;

    % Initial guesses
    rho_c = 1215; % initial crust density in kg/m^3
    rho_m = 2429; % initial mantle density in kg/m^3
    t_cr = 37.7e3; % initial reference crust thickness in meters

    % Airy isostasy
    dr = root_airy(rho_c, rho_m, h_topo);

    % perform transformation
    cs = GSHA(dr, 18);
    sc = cs2sc(cs);

    n_sh = 1:size(sc,1);

    phi = flexural_response(n_sh, g_avg, E_cr, T_e, sigma_cr, r_ceres, rho_c, rho_m);

    sc_flex = zeros(size(sc));

    for m = 1:size(sc,2)
        sc_flex(:,m) = sc(:,m).*phi';
    end

    dr_flex = GSHS(sc_flex,lonGrid,90-latGrid,18);

    % save intial data
    t_boundary_flex = - uniform_matrix .* t_cr - dr_flex;
    V_i = segment_2layer_model(h_topo, t_boundary_flex, -175e3, rho_c, rho_m, 20e3, Model);

    [n, DV] = degreeVariance(V_i);
    [n_obs, DV_obs] = degreeVariance(V_obs);

    res = sum((DV_obs(3:18) - DV(3:18)).^2 , "all");
end