function locmax = simulateLocMax(D, rho, var, niters, nondiag, nsubj)
% simulateLocMax(D, rho, var, niters) simulate the local maxima in an
% isotropic field through theoretical distribution of multivariate
% gaussian distribution.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the isotropic field.
% rho           spatial correlation between two adjacent voxels.
% var           variance of the field.
% niters        iteration times to generate the local maxima.
% nsubj         subject number used in multivariate t-statistics
% nondiag       whether we do not count the diagonal as the neighbor
%--------------------------------------------------------------------------
% OUTPUT
% the local maxima from theoretical multivariate gaussian distribution.
%--------------------------------------------------------------------------
% EXAMPLES
% D = 1;
% var = 1;
% rho = 0.8;
% niters = 10000;
% simulateLocMax(D, rho, var, niters)
%
if nargin == 4
    nsubj = 1000;
    nondiag = 0;
end

if nargin == 5
    nsubj = 1000;
end

locmax = zeros(1,niters);
k = 1;  % update the count of local maxima.
switch D
    case 1
        corr12 = rho;
        corr23 = rho^4;
        corr11 = 1;
        corr = [corr11, corr12, corr23; corr12,corr11, corr12;corr23, corr12, corr11];
        sigma = var*corr;
        dimMat = 3;
        midpt = 2;
    case 2
        corr12 = rho;
        corr23 = rho^4;
        corr11 = 1;
        corr1 = [corr11, corr12, corr23; corr12,corr11, corr12;corr23, corr12, corr11];
        corr = kron(corr1, corr1);
        sigma = var*corr;
        dimMat = 9;
        midpt = 5;
        diagind = sort([midpt-3.^(0:(D-1)), midpt+3.^(0:(D-1)), midpt]);
    case 3
        corr12 = rho;
        corr23 = rho^4;
        corr11 = 1;
        corr1 = [corr11, corr12, corr23; corr12,corr11, corr12;corr23, corr12, corr11];
        corr = kron(kron(corr1, corr1),corr1);
        sigma = var*corr;
        dimMat = 27;
        midpt = 14;
        diagind = sort([midpt-3.^(0:(D-1)), midpt+3.^(0:(D-1)), midpt]);
end

% using multivariate gaussian distribution
R = mvnrnd(zeros(1,dimMat), sigma, niters);
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


locmax = locmax(1:(k-1));

% using t-statistics
%{
R = mvnrnd(zeros(1,dimMat), Sigma, niters*nsubj);
for iter = 1:niters
    Rtstat = mvtstat(R((nsubj*(iter-1) + 1): nsubj*iter, :)');
    if Rtstat(midpt) == max(Rtstat)
        locmax(k) = max(Rtstat);
        k = k+1;
    end
end
locmax = locmax(1:(k-1));
%}
end