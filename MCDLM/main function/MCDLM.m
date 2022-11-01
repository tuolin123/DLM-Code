function [peaklocs, peakinds, peakvals, peakpvals] = MCDLM(lat_data,...
    mask, D, FWHM, isotropic, gaussianize, sigma, niter, nlim, nj, nondiag)
% MCDLM(lat_data, mask, D, FWHM, sigma, niter, nlim, nj, nondiag) find the
% top local maxima of an image that lies in a mask, and calculate the
% p-values of these local maxima by using MCDLM approach with 
% (1) t-distribution and (2) gaussian distribution after gaussianization.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data  a Dim by nsubj array of data
% mask      a 0/1 mask
% D         dimension of the image, length of Dim
% gaussianize  whether to apply gaussianization for the peak height, 1 for
%           apply and 0 for not, default is 0   
% FWHM      the FWHM of smoothing kernel applied to smooth the original
%           images
% isotropic a 0/1 value denotes for whether the user want to use the
%           isotropic property of the field, default is 1
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
%                   calculated for peakvals or gaussianized peakvals
%--------------------------------------------------------------------------
% EXAMPLES
%2D example (isotropic)
% Dim = [91,109];
% D = length(Dim);
% nimgs = 20;
% FWHM = 1.4;
% gaussianize = 0;
% isotropic = 1;
% lat_data = normrnd(0,1,[Dim,nimgs]);
% mask = zeros(Dim); mask(30:60,40:80) = 1; mask = logical(mask);
% [peaklocs, peakinds, peakvals, peakpvals] =...
% MCDLM(lat_data, mask, D, FWHM, isotropic, gaussianize);
%
%2D example (non isotropic)
% Dim = [91,109];
% D = length(Dim);
% nimgs = 20;
% FWHM = 1.4;
% gaussianize = 0;
% isotropic = 0;
% lat_data = normrnd(0,1,[Dim,nimgs]);
% mask = zeros(Dim); mask(30:60,40:80) = 1; mask = logical(mask);
% [peaklocs, peakinds, peakvals, peakpvals] =...
% MCDLM(lat_data, mask, D, FWHM, isotropic, gaussianize);
%
% %3D example (with gaussianization)
% Dim = [91,109,91];
% D = length(Dim);
% nimgs = 20;
% FWHM = 1.4;
% gaussianize = 1;
% isotropic = 1;
% lat_data = normrnd(0,1,[Dim,nimgs]);
% mask = zeros(Dim); mask(30:60,40:80,30:60) = 1; mask = logical(mask);
% [peaklocs, peakinds, peakvals, peakpvals] =...
% MCDLM(lat_data, mask, D, FWHM, isotropic, gaussianize);
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
if ~exist('FWHM', 'var')
    FWHM = 0;
end
if ~exist('isotropic', 'var')
    isotropic = 1;
end
if ~exist('gaussian', 'var')
    gaussianize = 0;
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
Dim = size(lat_data);
nimgs = Dim(D+1);
variable_index = repmat( {':'}, 1, D );
f_new = zeros(size(lat_data));

if FWHM == 0
    f_new = lat_data;
else
    for nn = 1:nimgs
        resadd = 0; enlarge = 0;
        params = ConvFieldParams( repelem(FWHM, D), resadd, enlarge);
        lattice_smooth = convfield(lat_data(variable_index{:}, nn), params);
        f_new(variable_index{:}, nn) = lattice_smooth.field;
    end
end


% Calculate one-sample t-field
lattice_tfield = convfield_t(f_new, 0); % mvtstat(f_new) 
[peaklocs, peakinds, peakvals] = lmindices(lattice_tfield.field,'all',mask);

% apply the mask
for nn = 1:nimgs
     f_new(variable_index{:}, nn) = squeeze(f_new(variable_index{:},nn)).*zero2nan(mask);
end

% standardize
f_new = f_new./std(f_new, 0, D+1);

% estimate the sd of smoothing kernel
    [~, sigma_est] = estsd(f_new, mask, D, 0.5); 

if isotropic == 1
    % transfer into the spatial correlation
    rho = discrete_covariance(mean(sigma_est), D);
    
    % run MCDLM with calculated rho
    if gaussianize == 0
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
    else
        % MCDLM method by applying the gaussianization
        locmaxZ2 = [];
        for j = 1:nj
            locmaxZ_sigma2 = simulateLocMax_Discrete(D,rho,sigma,niter,nondiag); %apply MCDLM to Gaussian field
            locmaxZ2 = [locmaxZ_sigma2 locmaxZ2];
            if length(locmaxZ2) > nlim
                break;
            end
        end
        
        peakvals_gauss = norminv(tcdf(peakvals, nimgs));
        [empf2, empz2] = ecdf(locmaxZ2);
        empz2 = empz2(2:end); empf2 = empf2(2:end);   %remove non-unique point
        peakpvals = 1-interp1(empz2,empf2,peakvals_gauss);
    end
else
    % transfer into the correlation
    rho = discrete_covariance_ellipsoid(sigma_est, D);
    
    % run MCDLM with calculated rho
    if gaussianize == 0
        locmaxZ1 = [];
        
        for j = 1:nj
            locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(D,rho,niter,nondiag,nimgs); %apply MCDLM to t-stats
            locmaxZ1 = [locmaxZ_sigma1 locmaxZ1];
            if length(locmaxZ1) > nlim
                break;
            end
        end
        
        [empf, empz] = ecdf(locmaxZ1);
        empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
        peakpvals = 1-interp1(empz,empf,peakvals);
    else
        % MCDLM method by applying the gaussianization
        locmaxZ2 = [];
        for j = 1:nj
            locmaxZ_sigma2 = simulateLocMax_Discrete_Stationary(D,rho,niter,nondiag,nimgs); %apply MCDLM to Gaussian field
            locmaxZ2 = [locmaxZ_sigma2 locmaxZ2];
            if length(locmaxZ2) > nlim
                break;
            end
        end
        
        peakvals_gauss = norminv(tcdf(peakvals, nimgs));
        [empf2, empz2] = ecdf(locmaxZ2);
        empz2 = empz2(2:end); empf2 = empf2(2:end);   %remove non-unique point
        peakpvals = 1-interp1(empz2,empf2,peakvals_gauss);
    end
end


end
