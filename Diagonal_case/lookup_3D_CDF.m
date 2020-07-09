%% Look-up table for 3D CDF value with different rho and x
%--------------------------------------------------------------------------
% rho: spatial correlation, varying from 0.01 to 0.99.
% x: value of local maxima 
% D: the dimension of the field, 3 in this look-up table.
% sigma: variance of the field.
% nx: number of local maximum we want.
% niter: the number of iterations in function simulateLocMax to get the
% number of local maximum we desire
% matcdf: the matrix of the look-up table, rho for row and x for column
drho = 0.01;
rho = drho:drho:(1-drho);
D = 3;
sigma = 1;
nx = 1e5;
max_iter = 1e7;
niter = round(nx./(4.6501e5-4.7875e5.*rho).*max_iter); % coefficient from linear regression of rho and sample iteration numbers
niter(96) = nx/7816*max_iter;
niter(97) = nx/5106*max_iter;
niter(98) = nx/2741*max_iter;
niter(99) = nx/1019*max_iter;
matcdf = zeros(1/drho-1, nx);

tic
for iter = 1:length(rho)
    % loop through all rho to get the local maxima values and its
    % corresponding empirical cdf
    locmaxZ = [];
    niter_valid = min(niter(iter),max_iter);
    locmaxZ1 = simulateLocMax(D,rho(iter),sigma,niter_valid);
    locmaxZ = [locmaxZ, locmaxZ1];
    [fx,x] = ecdf(locmaxZ); 
    x = x(2:end); fx = fx(2:end);   %remove non-unique point
    xarray{iter} = x; % local maxima value
    fxarray{iter} = fx; % empirical CDF value
end
toc
% 21.5 minutes to run
count3D = cellfun(@(x) numel(x),fxarray); 

% Sampling from the union for all x to get the location to interpolate
xarrayall = cat(1, xarray{:});
xlocmax = sort([randsample(xarrayall, nx-2);min(xarrayall);max(xarrayall)]);

% Interpolate for all rho
for iter = 1:length(rho)
    x = xarray{iter};
    fx = fxarray{iter};
    fxLocmax = interp1(x,fx,xlocmax);
    matcdf(iter,:) = fxLocmax;
end

% Replace all NaNs with 0 in the beginning and with 1 in the end
imp_begin = matcdf(:,1:1000);
TF_begin = isnan(imp_begin);
imp_begin(TF_begin) = 0;
matcdf(:,1:1000) = imp_begin;
imp_end = matcdf(:,(end-500):end);
TF_end = isnan(imp_end);
imp_end(TF_end) = 1;
matcdf(:,(end-500):end) = imp_end;
save('presmooth_lookup_3D.mat', 'matcdf');

%% smoothing the look-up table

% Polynomial regression (with more degrees) smoothing across rho
matcdf_poly_3D = matcdf;

% Cross-validation to find the best parameter of smoothing spline
lambda = 0.99:0.001:0.999;
lrho = length(rho);
seed = 2020;
rng(seed)
randInd = randsample(1:lrho, lrho);
randInd = sort(randInd);
randrho = rho(randInd);
randmat = matcdf(randInd,:);
kfold = 5;
startInd = 1:ceil(lrho/5):lrho;
mse = zeros(1,length(lambda));
randnx = randsample(1:nx, nx/10); % do CV only for nx/10 of them to save time

tic
for j = 1:length(lambda)
    %for i = randnx
    i=randnx(1);
        for k = 1:kfold % 5-fold CV
            if k ~= kfold
                X_test = randrho(startInd(k):startInd(k+1)-1)';
                y_test = randmat(startInd(k):startInd(k+1)-1,i);
                X_valid = randrho(setdiff(1:lrho, startInd(k):startInd(k+1)-1))';
                y_valid = randmat(setdiff(1:lrho, startInd(k):startInd(k+1)-1),i);
            else
                X_test = randrho(startInd(k):end)';
                y_test = randmat(startInd(k):end,i);
                X_valid = randrho(setdiff(1:lrho, startInd(k):lrho))';
                y_valid = randmat(setdiff(1:lrho, startInd(k):lrho),i);
            end
            yfit = csaps(X_valid,y_valid,lambda(j),X_test);
            mse(j) = mean((yfit-y_test).^2)+mse(j);
        end
    %end
    mse(j) = mse(j)/((nx/10)*kfold);
end
toc
[m1,I1] = min(mse);

tic
for i = 1:nx
    X = rho';
    y = matcdf(:,i);
    %pfit = polyfit(X, y, 3);
    %fitval = polyval(pfit,X);
    %pfit = fitlm(X, y, 'y~x1:x1:x1:x1:x1-1');
    %fitval = predict(pfit,X);
    fitval = csaps(X,y,lambda(I1),X); % lambda = 0.999
    matcdf_poly_3D(:,i) = fitval;
end
toc

% 1-D convolution across x
%conv_length = 31; % convolution interval length, needed odd
%pad = (conv_length-1)/2;
%convmat = ones(1,conv_length)*1/conv_length;
%matcdf_pad = [repmat(matcdf_poly(:,1),1,pad), matcdf_poly, repmat(matcdf_poly(:,end),1,pad)];
%matcdf_conv = conv2(matcdf_pad, convmat, 'valid');
% cubic spline smoothing
lx = length(xlocmax);
seed = 2020;
rng(seed)
randInd = randsample(1:lx, lx);
randInd = sort(randInd);
randx = xlocmax(randInd);
randmat = matcdf_poly_3D(:,randInd);
kfold = 5;
startInd = 1:ceil(lx/5):lx;
mse2 = zeros(1,length(lambda));
randrho = randsample(1:lrho, lrho/3);

tic
for j = 1:length(lambda)
    for i = 1:randrho
        %i=1;
        for k = 1:kfold % 5-fold CV
            if k ~= kfold
                X_test = randx(startInd(k):startInd(k+1)-1)';
                y_test = randmat(i,startInd(k):startInd(k+1)-1);
                X_valid = randx(setdiff(1:lx, startInd(k):startInd(k+1)-1))';
                y_valid = randmat(i,setdiff(1:lx, startInd(k):startInd(k+1)-1));
            else
                X_test = randx(startInd(k):end)';
                y_test = randmat(i,startInd(k):end);
                X_valid = randx(setdiff(1:lx, startInd(k):lx))';
                y_valid = randmat(i,setdiff(1:lx, startInd(k):lx));
            end
            yfit = csaps(X_valid,y_valid,lambda(j),X_test);
            mse2(j) = mean((yfit-y_test).^2)+mse2(j);
        end
    end
    mse2(j) = mse2(j)/((lrho/3)*kfold);
end
toc
[m2,I2] = min(mse2);

matcdf_conv_3D = matcdf;
tic
for i = 1:lrho
    X = xlocmax';
    y = matcdf_poly_3D(i,:);
    fitval = csaps(X,y,lambda(I2),X); % lambda = 0.99
    matcdf_conv_3D(i,:) = fitval;
end
toc

% save the final smoothed table
matcdf_smooth = matcdf_conv_3D;
cdf_table3D.p = matcdf_smooth;
cdf_table3D.rho = rho;
cdf_table3D.x = xlocmax;

save('lookup_3D_CDF.mat','cdf_table3D')

%% Why smoothing good
rdsample = randsample(1:100000, 50);
plot(rho,matcdf(:,rdsample),'LineWidth', 2)
xlabel('$\rho$','Interpreter','latex')
ylabel('$F(u;\rho)$','Interpreter','latex')
title('Pre-smooth')
set(gca,'FontSize', 18)

plot(rho,matcdf_smooth(:,rdsample),'LineWidth', 2)
xlabel('$\rho$','Interpreter','latex')
ylabel('$F(u;\rho)$','Interpreter','latex')
title('After-smooth')
set(gca,'FontSize', 18)

plot(rho, matcdf(:,10), rho, matcdf_smooth(:,10),'LineWidth', 2)
xlabel('$\rho$','Interpreter','latex')
ylabel('$F(u;\rho)$','Interpreter','latex')
title(['u = ' num2str(round(cdf_table3D.x(10),3))], 'Interpreter','latex')
set(gca,'FontSize', 21)
plot(rho, matcdf(:,50000), rho, matcdf_smooth(:,50000),'LineWidth', 2)
xlabel('$\rho$','Interpreter','latex')
ylabel('$F(u;\rho)$','Interpreter','latex')
title(['u = ' num2str(round(cdf_table3D.x(50000),3))], 'Interpreter','latex')
set(gca,'FontSize', 21)
plot(rho, matcdf(:,99990), rho, matcdf_smooth(:,99990),'LineWidth', 2)
xlabel('$\rho$','Interpreter','latex')
ylabel('$F(u;\rho)$','Interpreter','latex')
title(['u = ' num2str(round(cdf_table3D.x(99990),3))], 'Interpreter','latex')
set(gca,'FontSize', 21)

rdsample_rho = randsample(1:length(rho), 10);
plot(xlocmax,matcdf(rdsample_rho,:))
xlabel('u')
ylabel('$F(u;\rho)$','Interpreter','latex')
title('Pre-smooth')
set(gca,'FontSize', 18)

plot(xlocmax,matcdf_smooth(rdsample_rho,:))
xlabel('u')
ylim([0,1])
ylabel('$F(u;\rho)$','Interpreter','latex')
title('After-smooth')
set(gca,'FontSize', 18)

%% Plot of relationship between rho and number of local maxima
rho = exp(-1/2.*1./(2.*[0.6 1.3 2.1 3 3.8 5.1].^2));
numlocmax = [227483 45519 12423 4618 2198 979];
%scatter(rho, numlocmax);
fit = fitlm(rho, numlocmax);
plot(fit, 'LineWidth', 2);
xlabel('$\rho$','Interpreter','latex')
ylabel('# of local maxima')
title('3-D Case')
legend('# of local maxima after 1e7 iterations', 'linear fit', 'confidence bounds')
set(gca,'FontSize', 18)
