function [rho, sigma] = estsd(Z, mask, D, initial)
% estsd(Z, mask, D, initial) estimate the correlation of the adjacent 
% voxels in each of the lattice directions. Meanwhile, it can also returns
% the standard deviation of the smoothing filter.
%--------------------------------------------------------------------------
% ARGUMENTS
% Z        Data with size Dim, and with nan where there is missing
%          data
% mask     A mask of the data which is made of 1s and 0s. 1s for where
%          the data is inside the mask and 0s for where the data is ouside
%          the mask. Default is taken to be the mask with 1s everywhere
% D        Dimension of the data, varies from 1, 2, 3
% initial  The initial value for solving rho2sd() function
%--------------------------------------------------------------------------
% OUTPUT
% rho      The correlation of the adjacent voxels in each of the lattice
%          directions
% sigma    The standard deviation of the smoothing filter
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = [91,109];
% nsubj = 1;
% noise = noisegen(Dim, nsubj, 20);
% D = 2;
% initial = 0.5;
% mask = spm_read_vols(spm_vol('MNImask.nii'));
% mask_2d = squeeze(mask(:,:,45));
% estsd(noise, mask_2d, 2, 0.5)
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
Z = squeeze(Z);
szz = size(Z);
Dim = szz(1:D);

if ~exist('mask','var' )
    mask = ones(Dim(1:D));
end
if ~exist('D','var')
    D = length(Dim);
end
if ~isequal(size(mask), Dim)
    error('The mask must be the same size as the field')
end
if ~exist('initial','var' )
    initial = 0.5;
end

Z = Z.*zero2nan(mask);
Z_stand = Z/std(Z(~isnan(Z)));
Z_mean = mean(Z(~isnan(Z)));

switch D
    case 1
        rho_hor = (Z_stand(1:(Dim(1)-1),:)-Z_mean).*(Z_stand(2:Dim(1),:)-Z_mean);
        rho = max(mean(rho_hor(~isnan(rho_hor))));
        rho = min(0.99, rho);
        fun = @(sigma) rho2sd(rho, sigma, D);
        sigma = fzero(fun, initial);
    case 2
        rho_hor = (Z_stand(1:(Dim(1)-1),1:Dim(2),:)-Z_mean).*(Z_stand(2:Dim(1),1:Dim(2),:)-Z_mean);
        rho_vert = (Z_stand(1:Dim(1),1:(Dim(2)-1),:)-Z_mean).*(Z_stand(1:Dim(1),2:Dim(2),:)-Z_mean);
        rho = [max(mean(rho_hor(~isnan(rho_hor))),0) max(mean(rho_vert(~isnan(rho_vert))),0)];
        rho = min(0.99, rho);
        fun = @(sigma) rho2sd(rho(1), sigma, D);
        sigma1 = fzero(fun, initial);
        fun = @(sigma) rho2sd(rho(2), sigma, D);
        sigma2 = fzero(fun, initial);
        sigma = [sigma1 sigma2];
    case 3
        rho_hor = (Z_stand(1:(Dim(1)-1),1:Dim(2),1:Dim(3),:)-Z_mean).*(Z_stand(2:Dim(1),1:Dim(2),1:Dim(3),:)-Z_mean);
        rho_vert = (Z_stand(1:Dim(1),1:(Dim(2)-1),1:Dim(3),:)-Z_mean).*(Z_stand(1:Dim(1),2:Dim(2),1:Dim(3),:)-Z_mean);
        rho_depth = (Z_stand(1:Dim(1),1:Dim(2),1:(Dim(3)-1),:)-Z_mean).*(Z_stand(1:Dim(1),1:Dim(2),2:Dim(3),:)-Z_mean);
        rho = [max(mean(rho_hor(~isnan(rho_hor))),0) max(mean(rho_vert(~isnan(rho_vert))),0) max(mean(rho_depth(~isnan(rho_depth))),0)];
        rho = min(0.99, rho);
        fun = @(sigma) rho2sd(rho(1), sigma, D);
        sigma1 = fzero(fun, initial);
        fun = @(sigma) rho2sd(rho(2), sigma, D);
        sigma2 = fzero(fun, initial);
        fun = @(sigma) rho2sd(rho(3), sigma, D);
        sigma3 = fzero(fun, initial);
        sigma = [sigma1 sigma2 sigma3];
end