function [peaklocs, peakinds, peakvals, peakpvals, peakpvals_gauss] = MCDLM(lat_data,...
    mask, D, FWHM, sigma, niter, nlim, nj, nondiag)
% MCDLM(lat_data, mask, D, FWHM, sigma, niter, nlim, nj, nondiag) find the
% top local maxima of an image that lies in a mask, and calculate the
% p-values of these local maxima by using MCDLM approach with 
% (1) t-distribution and (2) gaussian distribution after gaussianization.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data  a Dim by nsubj array of data
% mask      a 0/1 mask
% D         dimension of the image, length of Dim
% FWHM      the FWHM of smoothing kernel applied to smooth the original
%           images
% sigma     the variance of the field, default is 1
% niter     the iteration time used for simulateLocMax_Discrete()
% nlim      the number of local maxima for the true height distribution
% nj        a number of simulation time to ensure the number of local
%           maxima generated for the true height distribution is enough,
%           this parameter can solve the memory problem when setting niter
%           low and this one high
% nondiag   a 0/1 value which denotes for whether we only consider the
%           partial connectivity case, default is 0
%--------------------------------------------------------------------------
% OUTPUT
% peaklocs          an D by npeaks array of the peak locations
% peakinds          an npeaks by 1 array where each value is the converted
%                   peak location, converted according to the size of 
%                   lat_data using convind
% peakvals          an npeaks by 1 array each value of which is the value
%                   the field takes at a given peak
% peakpvals         an npeaks by 1 array each value of which is the p-value
%                   calculated for peakvals
% peakpvals_gauss   an npeaks by 1 array each value of which is the p-value
%                   calculated for gaussianized peakvals
%--------------------------------------------------------------------------
% EXAMPLES
% %2D example
% Dim = [91,109];
% nimgs = 20;
% FWHM = 1.4;
% lat_data = normrnd(0,1,[Dim,nimgs]);
% mask = zeros(Dim); mask(30:60,40:80) = 1; mask = logical(mask);
% [peaklocs, peakinds, peakvals, peakpvals, peakpvals_gauss] =
% MCDLM(lat_data, mask, 2, FWHM);
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
if ~exist('FWHM', 'var')
    FWHM = 0;
end
if ~exist('sigma', 'var')
    sigma = 1;
end
if ~exist('niter', 'var')
    niter = 1e5;
end
if ~exist('nlim', 'var')
    nlim = 1e6;
end
if ~exist('nj', 'var')
    nj = 1e5; %increase the size of nj to solve memory issue
end
if ~exist('nondiag', 'var')
    nondiag = 0; %consider the diagonal
end

% Smoothing after reading the real data
if FWHM == 0
    f_new = lat_data;
else
    lattice_smooth = convfield(lat_data, FWHM);
    f_new = lattice_smooth.field;
end


% Calculate one-sample t-field
lattice_tfield = convfield_t(f_new, 0); 
[peaklocs, peakinds, peakvals] = lmindices(lattice_tfield.field,'all',mask);

% apply the mask
variable_index = repmat( {':'}, 1, D );
Dim = size(lat_data);
nimgs = Dim(D+1);

for nn = 1:nimgs
     f_new(variable_index{:}, nn) = squeeze(f_new(:,:,nn)).*zero2nan(mask);
end

% standardize
f_new = f_new./std(f_new, 0, D+1);

% estimate the correlation
[rho_est, sigma_est] = estsd(f_new, mask, D, 0.5); %if we assume stationarity
rho = discrete_covariance(mean(sigma_est), D);

% run MCDLM with calculated rho
locmaxZ1 = [];

for j = 1:nj
    locmaxZ_sigma1 = simulateLocMax_Discrete(D,rho,sigma,niter,nondiag,nimgs); %apply MCDLM to t-stats
    locmaxZ1 = [locmaxZ_sigma1 locmaxZ1];
    if length(locmaxZ1) > nlim
        break;
    end
end

[empf, empz] = ecdf(locmaxZ1);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
peakpvals = 1-interp1(empz,empf,peakvals);

% MCDLM method by applying the gaussianization
locmaxZ2 = [];

for j = 1:nj
    locmaxZ_sigma2 = simulateLocMax_Discrete(D,rho,sigma,niter); %apply MCDLM to Gaussian field
    locmaxZ2 = [locmaxZ_sigma2 locmaxZ2];
    if length(locmaxZ2) > nlim
        break;
    end
end

peakvals_gauss = norminv(tcdf(peakvals, nimgs));
[empf2, empz2] = ecdf(locmaxZ2);
empz2 = empz2(2:end); empf2 = empf2(2:end);   %remove non-unique point
peakpvals_gauss = 1-interp1(empz2,empf2,peakvals_gauss);

end

