function locmax = simulateLocMax(D, rho, var, niters, nondiag, nsubj)
% simulateLocMax(D, rho, var, niters) generates the local maxima in an
% isotropic field through simulating the theoretical distribution of
% multivariate gaussian distribution or t-distribution, using the
% covariance matrix calculated from a continuous field.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the isotropic field.
% rho           the spatial correlation between two adjacent voxels in a 
%               continuous random field
% var           the variance of the field
% niters        the number of iterations to generate the local maxima, note 
%               that the number of local maxima generated will be
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
% %3D fully connected example
% D = 3;
% var = 1;
% rho = 0.9;
% niters = 10000;
% simulateLocMax(D, rho, var, niters)
% 
% %3D partially connected example
% D = 3;
% var = 1;
% nu = 1.3;
% rho = 0.99;
% niters = 10000;
% nondiag = 1;
% simulateLocMax(D, rho, var, niters, nondiag)
%
% %2D t-field example
% D = 2;
% var = 1;
% nu = 1.3;
% rho = 0.9;
% niters = 10000;
% nondiag = 1;
% nsubj = 20;
% simulateLocMax(D, rho, var, niters, nondiag, nsubj)
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
if nargin == 4
    I_tstat = 0; % whether we use t-statistics
    nondiag = 0;
end

if nargin == 5
    I_tstat = 0;
    nsubj = 1000;
end

if nargin == 6
    I_tstat = 1;
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
        dimMat = 3; % number of voxels in the neighboring matrix
        midpt = 2; % index of midpoint of the vector
    case 2
        corr12 = rho;
        corr23 = rho^4;
        corr11 = 1;
        corr1 = [corr11, corr12, corr23; corr12,corr11, corr12;corr23, corr12, corr11];
        % the 2D correlation matrix is a kronecker product of two 1D correlation
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
        % the 3D correlation matrix is a kronecker product of three 1D correlation
        corr = kron(kron(corr1, corr1),corr1);
        sigma = var*corr;
        dimMat = 27;
        midpt = 14;
        diagind = sort([midpt-3.^(0:(D-1)), midpt+3.^(0:(D-1)), midpt]);
end

if I_tstat == 0
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
else
    % using t distribution
    R = (mvnrnd(zeros(1,dimMat), sigma, niters*nsubj))';
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