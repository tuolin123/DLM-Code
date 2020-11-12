function locmax = simulateLocMax_Discrete_Stationary(D, sigma2, niters, nondiag)
% simulateLocMax_Discrete_Stationary(D, rho, var, niters) simulate the local 
% maxima in an stationary gaussian field through theoretical distribution of
% multivariate gaussian distribution.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the stationary field, at least 2.
% cov           empirical covariance matrix estimated from the stationary 
%               field, returned from empiricalCov().
% niters        iteration times to generate the local maxima.
% nsubj         subject number used in multivariate t-statistics
% nondiag       whether we do not count the diagonal as the neighbor
%--------------------------------------------------------------------------
% OUTPUT
% the local maxima from theoretical multivariate gaussian distribution.
%--------------------------------------------------------------------------
%
if nargin == 3
    nondiag = 0;
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

% using multivariate gaussian distribution
R = mvnrnd(zeros(1,dimMat), sigma2, niters);
if nondiag == 0
    for iter = 1:niters
        if max(R(iter,:)) == R(iter,midpt) % R(iter,midpt) > 
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


locmax = locmax(1:(k-1));

end