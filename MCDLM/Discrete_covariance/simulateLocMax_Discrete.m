function locmax = simulateLocMax_Discrete(D, rho, var, niters, nondiag, nsubj)
% simulateLocMax(D, rho, var, niters) generates the local maxima in an
% isotropic field through simulating the theoretical distribution of
% multivariate gaussian distribution or t-distribution, using the
% covariance matrix calculated from the field in a discrete lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the isotropic field
% rho           a vector of spatial correlation between all pairs of 
%               neighboring voxels with distinct values when the field is 
%               in a discrete lattice, generated from discrete_covariance.m
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
% nu = 1.3;
% rho = discrete_covariance(nu, D);
% niters = 10000;
% simulateLocMax_Discrete(D, rho, var, niters)
% 
% %3D partially connected example
% D = 3;
% var = 1;
% nu = 1.3;
% rho = discrete_covariance(nu, D);
% niters = 10000;
% nondiag = 1;
% simulateLocMax_Discrete(D, rho, var, niters, nondiag)
%
% %2D t-field example
% D = 2;
% var = 1;
% nu = 1.3;
% rho = discrete_covariance(nu, D);
% niters = 10000;
% nondiag = 1;
% nsubj = 20;
% simulateLocMax_Discrete(D, rho, var, niters, nondiag, nsubj)
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
        I = 1:3;
        f = @(x,y)(x-y).^2;
        euc = bsxfun(f,I,I'); % calculating the square of euclidean distance
        corr = zeros(3,3);
        corr(euc == 0) = rho(1); % assignin the rho values when square of euclidean distance is 0,1 or 4
        corr(euc == 1) = rho(2);
        corr(euc == 4) = rho(3);
        sigma = var*corr;
        dimMat = 3; % number of voxels in the neighboring matrix
        midpt = 2; % index of midpoint of the vector
    case 2
        sz = [3 3];
        [row, col] = ind2sub(sz,1:9);
        f = @(x,y)(x-y).^2;
        euc = bsxfun(f,row,row') + bsxfun(f,col,col');
        corr = zeros(9,9);
        corr(euc == 0) = rho(1);
        corr(euc == 1) = rho(2);
        corr(euc == 2) = rho(3);
        corr(euc == 4) = rho(4);
        corr(euc == 5) = rho(5);
        corr(euc == 8) = rho(6);
        sigma = var*corr;
        dimMat = 9;
        midpt = 5;
        diagind = sort([midpt-3.^(0:(D-1)), midpt+3.^(0:(D-1)), midpt]); 
        % the location index of the diagonals in a neighboring matrix
    case 3
        sz = [3 3 3];
        [I1, I2, I3] = ind2sub(sz, 1:27);
        f = @(x,y)(x-y).^2;
        euc = bsxfun(f,I1,I1') + bsxfun(f,I2,I2') + bsxfun(f,I3,I3');
        corr = zeros(27,27);
        corr(euc == 0) = rho(1);
        corr(euc == 1) = rho(2);
        corr(euc == 2) = rho(3);
        corr(euc == 3) = rho(4);
        corr(euc == 4) = rho(5);
        corr(euc == 5) = rho(6);
        corr(euc == 6) = rho(7);
        corr(euc == 8) = rho(8);
        corr(euc == 9) = rho(9);
        corr(euc == 12) = rho(10);
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