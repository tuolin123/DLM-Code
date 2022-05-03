function [rho, sigma] = estsd(Z, mask, D, initial)
% estsd(Z, mask, D, initial) estimates the correlation of the adjacent 
% voxels in each of the lattice directions. Meanwhile, it can also return
% the standard deviation of the smoothing filter.
%--------------------------------------------------------------------------
% ARGUMENTS
% Z        Data with size Dim * n, where Dim is the image size and n is the
%          sample size of images. Missing data is denoted as nan
% mask     A mask (size Dim) of the data which is made of 1s and 0s. 1s for where
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
% nsubj = 100;
% FWHM = 2;
% noise = noisegen(Dim, nsubj, FWHM);
% D = 2;
% initial = 0.5;
% mask = spm_read_vols(spm_vol('MNImask.nii'));
% mask_2d = squeeze(mask(:,:,45));
% [rho_est, sigma_est] = estsd(noise, mask_2d, 2, initial);
% sigma2FWHM(sigma_est) % check FWHM
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
Z = squeeze(Z);
z_size = size(Z);
Dim = z_size(1:D);

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
    initial = 0.1;
end

if length(z_size) == D
    Z = Z.*zero2nan(mask);
    Z_stand = Z./std(Z(~isnan(Z)));
    Z_mean = mean(Z_stand(~isnan(Z_stand)));
else
    Z_stand = Z./std(Z,0,D+1);
    Z_stand = Z_stand.*zero2nan(mask);
    Z_mean = mean(Z_stand(~isnan(Z_stand)));
end

switch D
    case 1
        rho_hor = (Z_stand(1:(Dim(1)-1),:)-Z_mean).*(Z_stand(2:Dim(1),:)-Z_mean);
        rho = max(mean(rho_hor(~isnan(rho_hor))),0);
    case 2
        rho_hor = (Z_stand(1:(Dim(1)-1),1:Dim(2),:)-Z_mean).*(Z_stand(2:Dim(1),1:Dim(2),:)-Z_mean);
        rho_vert = (Z_stand(1:Dim(1),1:(Dim(2)-1),:)-Z_mean).*(Z_stand(1:Dim(1),2:Dim(2),:)-Z_mean);
        rho = [max(mean(rho_hor(~isnan(rho_hor))),0) max(mean(rho_vert(~isnan(rho_vert))),0)];
    case 3
        rho_hor = (Z_stand(1:(Dim(1)-1),1:Dim(2),1:Dim(3),:)-Z_mean).*(Z_stand(2:Dim(1),1:Dim(2),1:Dim(3),:)-Z_mean);
        rho_vert = (Z_stand(1:Dim(1),1:(Dim(2)-1),1:Dim(3),:)-Z_mean).*(Z_stand(1:Dim(1),2:Dim(2),1:Dim(3),:)-Z_mean);
        rho_depth = (Z_stand(1:Dim(1),1:Dim(2),1:(Dim(3)-1),:)-Z_mean).*(Z_stand(1:Dim(1),1:Dim(2),2:Dim(3),:)-Z_mean);
        rho = [max(mean(rho_hor(~isnan(rho_hor))),0) max(mean(rho_vert(~isnan(rho_vert))),0) max(mean(rho_depth(~isnan(rho_depth))),0)];
end

sigma = rho2sd(rho, D, initial);

end