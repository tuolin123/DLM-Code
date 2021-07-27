function g = gaussian_kernel(D, sigma, x)
% Computes density function of the isotropic Gaussian kernel.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the kernel
% sigma         the standard deviation of the kernel
% x             an array for which the density of the Gaussian function
%               should be evaluated
%--------------------------------------------------------------------------
% OUTPUT
% g             the density of Gaussian function evaluated at x
%--------------------------------------------------------------------------
% EXAMPLES
% D = 2
% sigma = 1.3
% x = -5:1:5
% gaussian_kernel(D, sigma, x)
%__________________________________________________________________________
g = 1/(sqrt(2*pi)*sigma)^D*exp(-x.^2./(2*sigma^2));
end