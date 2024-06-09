clear;
close all;
clc;

% parameters
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
rho_c = 1200; 
rho_m = 2000; 
degree = 179;

latLim = [-89.5 89.5 1];
lonLim = [-179.5 179.5 1];
lonGrid = lonLim(1):lonLim(3):lonLim(2);
latGrid = latLim(1):latLim(3):latLim(2);
Lon_matrix = repmat(lonGrid,length(latGrid),1);
Lat_matrix = repmat(latGrid',1,length(lonGrid));
[num_rows, num_cols] = size(Lat_matrix);
uniform_matrix = ones(size(Lat_matrix));

% Set up synthetic topographic profiles
amplitude = 2.5e3;
frequency = 0.25; 

[x, y] = meshgrid(1:num_cols, 1:num_rows);

% Create the surface using a combination of sine and cosine functions
h_topo = 0.5 * amplitude * (sin(frequency * y) + cos(frequency * x));

dr = root_airy(rho_c, rho_m, h_topo);

%% GSHA

cs = GSHA(dr, degree);
sc = cs2sc(cs);

n_sh = 1:size(sc,1);

phi = flexural_response(n_sh, g_avg, E_cr, T_e_cr, sigma_cr, r_ceres, rho_c, rho_m);

sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*phi';
end

%% GSHS

dr_flex = GSHS(sc_flex,lonGrid,90-latGrid,179);

dr_flex = Europe_centered(dr_flex);

%%

figure('Position',[100 100 800 400]);
imagesc(lonGrid, latGrid(1+gt:end-gt), dr(1+gt:end-gt,:) ./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Root (km)');
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
imagesc(lonGrid, latGrid(1+gt:end-gt), dr_flex(1+gt:end-gt,:) ./1e3);
colorbar;
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Root with flexure response (km)');
colormap(turbo); 
set(gca, 'YDir', 'normal');
set(gca, 'XTick', -180:45:180);
set(gca, 'YTick', -90:45:90);
set(gca, 'XTickLabel', -180:45:180);
set(gca, 'YTickLabel', -90:45:90);
set(gca, 'XMinorTick', 'on', 'XMinorGrid', 'off');
set(gca, 'YMinorTick', 'on', 'YMinorGrid', 'off');
set(gca, 'TickDir', 'out');

% define functions
function dr = root_airy(rho_c, rho_m, dh)

    dr = (rho_c /  (rho_m - rho_c)) .* dh;

end

function phi_n = flexural_response(n, g, E, T_e, sigma, R, rho_c, rho_m)
    
    D = (E * T_e^3) / (12 * (1 - sigma^2));

    phi_n = (1 + (D / (g .* (rho_m - rho_c))) .* (2 .* (n + 1) ./ (2*R)).^4).^(-1);

end