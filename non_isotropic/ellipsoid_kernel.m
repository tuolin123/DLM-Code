function e = ellipsoid_kernel(D, sigma, x)
% Computes density function of the elliptical Gaussian kernel.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the kernel
% sigma         the standard deviation of the kernel in each direction
% x             an array for which the density of the Gaussian function
%               should be evaluated
%--------------------------------------------------------------------------
% OUTPUT
% g             the density of the elliptical Gaussian function evaluated 
%               at x
%--------------------------------------------------------------------------
% EXAMPLES
% #2D example
% D = 2;
% sigma = [1.3 1.3];
% [X, Y] = meshgrid((-5:1:5) ,(-5:1:5));
% gaussian_kernel(D, sigma, {X,Y})
%__________________________________________________________________________
x2 = 0;        
for i = 1:D
    x2 = x2 + cell2mat(x(i)).^2./(2*sigma(i)^2);
end
e = 1/(sqrt(2*pi))^D*prod(sigma)*exp(-x2);
end