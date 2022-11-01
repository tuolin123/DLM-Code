function locmax = simulateLocMax_Discrete_Stationary(D, sigma2, niters, nondiag, nsubj)
% simulateLocMax_Discrete_Stationary(D, sigma2, niters, nondiag) generates 
% the local maxima in an stationary gaussian field through simulating the
% theoretical distribution of multivariate gaussian distribution or 
% t-distribution, using the covariance matrix either from the theoretical 
% or empirical calculation.
%
% This is a more general function than simulateLocMax.m and 
% simulateLocMax_Discrete.m in that it allows any types of covariance
% function, but the whole covariance matrix is needed.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the stationary field, at least 2.
% sigma2        the covariance matrix estimated from the stationary 
%               field, returned from empiricalCov() if it's estimated 
%               empirically.
% niters        the number of iterations to generate the local maxima,
%               note that the number of local maxima generated will be
%               smaller than this number since we throw away those not
%               qualified as local maxima
% nondiag       a 0/1 value which denotes for whether we only consider the
%               partial connectivity case, default is 0
% nsubj         the degrees of freedom used in multivariate t-statistics, 
%               if specified, the function will use t-distribution
%--------------------------------------------------------------------------
% OUTPUT
% locmax        a set of simulated local maxima
%--------------------------------------------------------------------------
% EXAMPLES
% %apply the method using empirical covariance function
% Dim = [50 50]; nsubj = 100; D = length(Dim); FWHM = 2;
% noise = noisegen(Dim, nsubj, FWHM);
% sigma2 = empiricalCov(D, noise);
% 
% % performs bad with small sample size, need to ensure covariance spd
% [V,U] = eig(sigma2);
% U = max(10e-10, U);
% sigma2 = V*U*V';
% 
% niters = 10000;
% nondiag = 1;
% s = simulateLocMax_Discrete_Stationary(D, sigma2, niters, nondiag);
% 
%
% %apply the method using covariance function by calculation
% nu = [0.6 1.3];
% D = 2;
% sigma2 = discrete_covariance_ellipsoid(nu, D);
% 
% niters = 10000;
% nondiag = 1;
% s = simulateLocMax_Discrete_Stationary(D, sigma2, niters, nondiag);
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------

if nargin == 3
    I_tstat = 0; % whether we use t-statistics
    nondiag = 0;
end

if nargin == 4
    I_tstat = 0;
end

if nargin == 5
    I_tstat = 1;
end

locmax = zeros(1,niters);
k = 1;  % update the count of local maxima.
switch D
    case 2
        dimMat = 9;
        midpt = 5;
        diagind = sort([midpt-3.^(0:(D-1)), midpt+3.^(0:(D-1)), midpt]);
    case 3
        dimMat = 27;
        midpt = 14;
        diagind = sort([midpt-3.^(0:(D-1)), midpt+3.^(0:(D-1)), midpt]);
end

if I_tstat == 0
    % using multivariate gaussian distribution
    R = mvnrnd(zeros(1,dimMat), sigma2, niters);
    if nondiag == 0
        for iter = 1:niters
            if max(R(iter,:)) == R(iter,midpt) 
                locmax(k) = R(iter,midpt);
                k = k + 1;
            end
        end
    elseif nondiag == 1
        for iter = 1:niters
            if max(R(iter,diagind)) == R(iter,midpt)
                locmax(k) = R(iter,midpt);
                k = k + 1;
            end
        end
    end
else
    % using t distribution
    R = (mvnrnd(zeros(1,dimMat), sigma2, niters*nsubj))';
    R = reshape(R,[dimMat,niters,nsubj]);
    Rtstat = mvtstat(R,[dimMat,niters]);
    for iter = 1:niters
        if Rtstat(midpt,iter) == max(Rtstat(:,iter))
            locmax(k) = max(Rtstat(:,iter));
            k = k+1;
        end
    end
end



locmax = locmax(1:(k-1));

end