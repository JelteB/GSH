function [U,L,R,A] = import_layer_NEW(Model,nlayer)
%
% Fixed cus the previous code was hot garbooooo

% name of the layer
layer_name = ['l' num2str(nlayer)];
layer_base = ['l' num2str(nlayer+1)];

% filenames
fupper = Model.(layer_name).bound;
flower = Model.(layer_base).bound;
fdens = Model.(layer_name).dens;

U = gmt2matrix(fupper).*1e3;
L = gmt2matrix(flower).*1e3;

% Density file
R = gmt2matrix(fdens).*1e3;

