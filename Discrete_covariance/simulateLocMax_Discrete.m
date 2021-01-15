function locmax = simulateLocMax_Discrete(D, rho, var, niters, nondiag, nsubj)
% simulateLocMax(D, rho, var, niters) simulate the local maxima in an
% isotropic field through theoretical distribution of multivariate
% gaussian distribution.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the isotropic field.
% rho           spatial correlation between two voxels, a vector with all
% pairs in the neighborhood generated from discrete_covariance.
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
% rho = discrete_covariance(1.3, D);
% niters = 10000;
% simulateLocMax(D, rho, var, niters)
%
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
        dimMat = 3;
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
        %R = (mvnrnd(zeros(1,dimMat), sigma, nsubj))';
        %Rtstat = mvtstat(R,[dimMat,1]);
        if Rtstat(midpt,iter) == max(Rtstat(:,iter))
        %if Rtstat(midpt) == max(Rtstat)
            locmax(k) = max(Rtstat(:,iter));
            %locmax(k) = max(Rtstat);
            k = k+1;
        end
    end
end


locmax = locmax(1:(k-1));

end