function sd = rho2sd(rho, D, initial)
% Compute the standard deviation of the smoothing kernel based on the
% correlation between two adjacent voxels for a isotropic Gaussian field on
% a discrete lattice. Note this is different from the one for continuous
% field!!!
%--------------------------------------------------------------------------
% ARGUMENTS
% rho           a vector of correlations between two adjacent voxels in 
%               each lattice direction
% D             the dimension of the field, default is 1
% initial       initial value used to calculate the root, default is 0.1
%--------------------------------------------------------------------------
% OUTPUT
% sd            the standard deviation of the smoothing kernel in each
%               lattice direction
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D example
% rho2sd(0.9)
%
% % 2D example
% rho = [0.9 0.1];
% D = 2;
% rho2sd(rho, D)
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
switch nargin
    case 1
        D = 1;
        initial = 0.1;
    case 2
        initial = 0.1;
end

switch D
    case 1
        fun = @(sigma) diff(rho, sigma, D);
        sd = fzero(fun, initial);
    case 2
        fun = @(sigma) diff(rho(1), sigma, D);
        sigma1 = fzero(fun, initial);
        fun = @(sigma) diff(rho(2), sigma, D);
        sigma2 = fzero(fun, initial);
        sd = [sigma1 sigma2];
    case 3
        fun = @(sigma) diff(rho(1), sigma, D);
        sigma1 = fzero(fun, initial);
        fun = @(sigma) diff(rho(2), sigma, D);
        sigma2 = fzero(fun, initial);
        fun = @(sigma) diff(rho(3), sigma, D);
        sigma3 = fzero(fun, initial);
        sd = [sigma1 sigma2 sigma3];
end

    % this is a nested function    
    function d = diff(rho, sigma, D)
        cov = discrete_covariance(sigma, D);
        d = rho - cov(2);
    end
end