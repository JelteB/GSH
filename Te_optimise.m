clear;
close all;
clc;
HOME = pwd;
format long
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% % Define the objective function to minimize the difference in degree variance
objective_function = @(T_e) compute_DV_difference(T_e);

% Find the optimal T_e that minimizes the difference in degree variance
optimal_Te = fminbnd(@compute_DV_difference, 5, 1e9);

% Display the optimal T_e
fprintf('Optimal value of Te: %.2f\n', optimal_Te);

% Function to compute the difference in degree variance for a given Te
function error = compute_DV_difference(T_e)
    % Update the flexural response with the current value of T_e
    % define parameters
    r_ceres = 470e3;
    sigma_cr = 0.5;
    E_cr = 5e9;
    T_e_cr = T_e;
    g_avg = 0.27;  % m/s²

    % load data
    [topo_data] = load("Topo_Ceres.mat");

    % process data
    h_topo_raw = topo_data.Topography;
    h_topo = zeros(size(h_topo_raw));
    h_topo(:,end/2+1:end) = h_topo_raw(:,1:end/2);
    h_topo(:,1:end/2) = h_topo_raw(:,end/2+1:end);

    % Initial guesses
    rho_c = 1215; % initial crust density in kg/m^3
    rho_m = 2429; % initial mantle density in kg/m^3

    % Airy isostasy
    dr = (rho_c /  (rho_m - rho_c)) .* h_topo;

    % perform transformation
    cs = GSHA(dr, 18);
    sc = cs2sc(cs);
    n_sh = 1:size(sc,1);
    phi = flexural_response(n_sh, g_avg, E_cr, T_e_cr, sigma_cr, r_ceres, rho_c, rho_m);
    sc_flex = zeros(size(sc));

    for m = 1:size(sc,2)
        sc_flex(:,m) = sc(:,m).*phi';
    end

    % observation SH coeffs
    obs_data = load("SHcoeffs_Ceres.mat");
    Obs = obs_data.SHcoeff;
    [~, observed_DV] = degreeVariance(Obs);

    % Compute the degree variance for the updated model SH coefficients
    [Clm,Slm,llvec,mmvec] = sc2vecml(sc_flex,18);
    V_flex = [llvec' mmvec' Clm Slm];
    [~, model_DV] = degreeVariance(V_flex);
        
    error = sum(abs(observed_DV(4:9) - model_DV(4:9)));
end
% function
function phi_n = flexural_response(n, g, E, T_e, sigma, R, rho_c, rho_m)
    
    D = (E * T_e^3) / (12 * (1 - sigma^2));

    phi_n = (1 + (D ./ (g .* (rho_m - rho_c))) .* (2 .* (n + 1) ./ (2*R)).^4).^(-1);

end
