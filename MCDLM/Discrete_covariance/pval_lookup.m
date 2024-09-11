function pval = pval_lookup(rho, z, D)
% pval_lookup(D, rho, z) computes the p_value of a perticular height using lookup table method by interpolation.
%--------------------------------------------------------------------------
% ARGUMENTS
% z             a vector of local maxima value.
% rho           spatial correlation, dimension 3, in each direction lattice. Input should be a n*d matrix with n voxels and d dimensions.
% D             the dimension of the image.
%--------------------------------------------------------------------------
% OUTPUT
% the p-value of the  at value z.
%--------------------------------------------------------------------------
switch D
    case 1
        load('presmooth_lookup_1D_discretecov.mat')
        load('lookup_1D_CDF_discretecov.mat')
    case 2
        load('presmooth_lookup_2D_discretecov.mat')
        load('lookup_2D_CDF_discretecov.mat')
    case 3
        load('presmooth_lookup_3D_discretecov.mat')
        load('lookup_3D_CDF_discretecov.mat')
end

dp = cdf_table2D.p;
dx = cdf_table2D.x;
drho = cdf_table2D.rho;
dp(dp<0) = 0;
[X,Y] = meshgrid(dx,drho);

pval = 1-interp2(X,Y,dp,z,rho);
return
