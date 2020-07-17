function g = gaussian_kernel(D, sigma, x)
%__________________________________________________________________________
% Computes values of the Gaussian kernel.
% 
% Input:
%   x    -    Array for which the values of the Gaussian function should be
%             evaluated
%   D    -    The dimension of the kernel
%   sigma-    The standard deviation of the Gaussian kernel
%
% Output:
%   y:   -    Array of same size of x containing the values of the Gaussian
%             kernel
%__________________________________________________________________________

g = 1/(sqrt(2*pi)*sigma)^D*exp(-x.^2./(2*sigma^2));
end