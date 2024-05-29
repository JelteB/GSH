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

% In case there is a linear density gradient in the model
if nargout == 4
    if isfield(Model.(layer_name), 'alpha')
        falpha = Model.(layer_name).alpha; 
        if ischar(falpha)
            A = gmt2matrix(load(falpha)).*1e3;
        elseif isnumeric(falpha)
            A = falpha;
            if isscalar(falpha)
                A = falpha.*ones(size(U));
            end
        else
            error(['Could not resolve linear density gradient class & size:  ' class(falpha)])
        end         
    else
        error('No specific file is listed for the linear gradient in the layer')
    end
end