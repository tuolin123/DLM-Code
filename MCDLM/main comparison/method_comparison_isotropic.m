%% This script compares the performance of MCDLM, ADLM and continuous 
% methods for computing the height distribution of local maxima in a 2D(3D)
% isotropic gaussian random field.

%% Parameter set-up
D = 2; % dimension of the random field
sigma = 1; % variance of the random field
niter = 1e5; % pre-specified iteration times

% use the following line if the covariance function is derived from the 
% discrete lattice
rho_disc = 0.01;
if rho_disc == 0.01
    n_lim = 1e6;
elseif rho_disc == 0.5
    n_lim = 1e6;
elseif rho_disc == 0.99
    n_lim = 2e5;
end
% use the following line if the covariance function is derived from the
% continuous random field, change rho_disc to rho_cont for all code below
% stddev_fwhm = 5;
% rho_cont = round(exp(-1/2*1/(2*stddev_fwhm^2)),2);

nondiag = 0; % whether voxels in diagonal direction are considered as neighbors
% use the following line if the voxels in the diagonal direction are not
% considered as neighbors
% nondiag = 1;



%% MCDLM method implementation
% calculate the standard deviation of the smoothing kernel corresponds to
% the correlation of the field; skip this if rho_cont used
load('inverse_rho_discrete.mat')
stddev_fwhm = invrho(rho_disc*100);

% we now have a function to calculate the fwhm from rho_disc concisely
% stddev_fwhm = rho2sd(rho_disc, D);

tic
locmaxZ1 = [];
rho = discrete_covariance(stddev_fwhm, D); 
for j = 1:10000 % for loop implemented to save memory in each loop
    % If diagonal voxels are not considered, nondiag = 1; otherwise nondiag
    % = 0
    locmaxZ = simulateLocMax_Discrete(D,rho,sigma,niter,nondiag);
    
    % use below if rho_cont is used
    % locmaxZ = simulateLocMax(D,rho,sigma,niter,nondiag);
    
    % use below for estimation of rho from field
    % cov = empiricalCov(D,field); % field is the data you want to estimate the spatial correlation for
    % Make the covariance matrix spd
    % [V,DD] = eig(cov);
    % DD = max(10e-10, DD);
    % sigma2 = V*DD*V';
    % locmaxZ = simulateLocMax_Discrete_Stationary(D,sigma2,niter,nondiag);
    % alternatively, if isotropy is assured, we can use
    % locmaxZ = simulateLocMax_Discrete(D,sigma2(2,1),sigma,niter,nondiag);
    
    locmaxZ1 = [locmaxZ locmaxZ1];
    if length(locmaxZ1) > n_lim
        break;
    end
end
toc

%% Continuous method: density of local maxima(continuous case)
kappa   = 1; % theoretical kappa
density = peakHeightDensity( D, kappa );
dz = 0.001;
z=-2.5:dz:5; 
qcont = density(z);
%% ADLM method: density of local maxima(discrete case)
nz=length(z);
ddrho = 0.001;
drho=0:ddrho:0.999;
nrho=length(drho);
p=zeros(nz,nrho);
nnb=ones(nz,D)*2;
q=ddlm(z,ones(nz,D)*rho_disc,nnb,D);
qnormalize=q./(sum(q)*dz);

%% Simulate the fields and obtain the local maxima, try to get the true 
% height distribution of local maxima in a certain field
locmaxZ2 = [];
nsim = 1000;
dim = [50 50];
f_new = zeros([dim,nsim]);
cut = 1;
FWHM = 2*sqrt(2*log(2))*stddev_fwhm;
edge = ceil(4*stddev_fwhm);

for j = 1:10
    for nn = 1:nsim
        lat_data = normrnd( 0, 1, dim(1)+2*edge, dim(2)+2*edge);
        smoothed_fconv = fconv( lat_data, FWHM, 2);
        f_new(:,:,nn) = smoothed_fconv((edge+1):(edge+dim(1)),(edge+1):(edge+dim(2)));
        %variable_index = repmat( {':'}, 1, D );
        %f_new(variable_index{:}, nn)
    end

    f_new = f_new/std(f_new(:));
    
    for nn = 1:nsim
        Z = f_new(:,:,nn);
        Z = squeeze(Z);
        
        % Find local maxima of the field Z and remove maxima at the boundary
        Imax = imregionalmax(Z); Imin = imregionalmin(Z);
        % If not including the diagonal as neighbors, use the next line
        % Imax = imregionalmax(Z,4); Imin = imregionalmin(Z,4);
        
        % Not including the boundary
        Imax = Imax((1+cut):(end-cut), (1+cut):(end-cut));
        Imin = Imin((1+cut):(end-cut), (1+cut):(end-cut));
        
        % Rescale Z to match the dimensions of Imax and Imin
        Z = Z((1+cut):(end-cut), (1+cut):(end-cut));
        
        % add the new local maxima of Z to the vector containing all local
        % maxima
        locmaxZ2 = [locmaxZ2; Z(Imax); -Z(Imin)]; 
    end
end

%% Draw the pp-plot for comparison
% If using the pre-saved lookup table, use the next line
% tic
% pval_dist = pval_lookup(rho_disc, z, D, 'discrete',1); 
% toc

% Calculate the p-value by MCDLM
[empf, empz] = ecdf(locmaxZ1);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z);

figure();
% Calculate the reference p-value 
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

% Plot the pp-plot
plot1 = plot(pval_dist,pval_ref,'LineWidth', 2);
plot1;
xlim([0 0.05])
ylim([0 0.05])
hold on
pval_cont = 1-cumsum(qcont).*dz;
plot2 = line(pval_cont,pval_ref, 'Color','r', 'LineWidth', 2);
plot2;
hold on
pval_disc = 1-cumsum(qnormalize).*dz;
plot3 = line(pval_disc,pval_ref, 'Color', 'green', 'LineWidth', 2);
plot3;
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
set(gca,'FontSize', 15)
xlabel('$p$-value', 'Interpreter', 'latex', 'fontsize',18)
ylabel('Reference $p$-value', 'Interpreter', 'latex', 'fontsize',18)
xticks(0:0.01:0.05)
%lgd = legend([plot2, plot1, hline], 'Continuous RFT','MCDLM', '45 degree line', 'Location','northwest');
lgd = legend([plot3, plot1, plot2, hline], 'ADLM', 'MCDLM', 'Continuous RFT','45 degree line', 'Location','southeast');
lgd.FontSize = 12;
title(['$\rho$ = ' num2str(rho_disc) '\ (' 'FWHM = ' num2str(round(FWHM,1)) ')'], 'Interpreter','latex', 'fontsize',18)
axis square

%saveas(gcf, 'lookup_2D_rho0.01_disc.jpg')
%% Plot using Armin and Fabian's method
[empf, empz] = ecdf(locmaxZ2);
empz = empz(2:end); empf = empf(2:end);
z = empz;
dz = diff(z);
pval_ref = 1-interp1(empz,empf,z); 
[empp, empzz] = ecdf(pval_ref);
empp = empp(2:end); empzz = empzz(2:end);
plot2 = plot(empzz,empp,'--','Color','k');
plot2;
hold on

[empf, empz] = ecdf(locmaxZ1);
empz = empz(2:end); empf = empf(2:end);
%pval_dist = 1-empf;
pval_dist = 1-interp1(empz,empf,z);
[empp, empzz] = ecdf(pval_dist);
empp = empp(2:end); empzz = empzz(2:end);
plot1 = plot(empzz,empp,'Color',[0, 0.4470, 0.7410],'LineWidth', 2);
plot1;
hold on

qcont = density(z);
pval_cont = 1-cumsum(qcont(2:end,:).*dz);
[empp, empzz] = ecdf(pval_cont);
empp = empp(2:end); empzz = empzz(2:end);
plot4 = plot(empzz,empp,'r','LineWidth', 2);
plot4;
hold on

% nz = length(z);
% nnb = ones(nz,3)*2;
% q=ddlm(z,ones(nz,3)*rho_disc,nnb,2);
% qnormalize=q./(sum(q)*dz);
% pval_disc = 1-cumsum(qnormalize).*dz;
% [empp, empzz] = ecdf(pval_disc);
% empp = empp(2:end); empzz = empzz(2:end);
% plot5 = plot(empzz,empp,'g','LineWidth', 2);
% plot5;
% hold on

xlim([0 0.05])
ylim([0 0.05])
set(gca,'FontSize', 18)
xlabel('$p$-value', 'Interpreter', 'latex', 'fontsize',24)
ylabel('Reference $p$-value', 'Interpreter', 'latex', 'fontsize',24)
legend([plot4,plot1,plot2], 'Continuous RFT','MCDLM', '45 degree line', 'Location','northwest')
title(['$\rho$ = ' num2str(rho_disc) '\ (' 'FWHM = ' num2str(round(FWHM,1)) ')'], 'Interpreter','latex', 'fontsize',24)
%ax = gca;
axis square

%% Using the empirical covariance function
% tic
% locmaxZ =[];
% for i = 1:1e3
% rho_vec = discrete_covariance(invrho(ceil(100*rho_disc)), D);
% locmaxZ_new = simulateLocMax_Discrete(D,rho_vec,sigma,1e4);
% locmaxZ = [locmaxZ, locmaxZ_new];
% end
% toc
sigma2 = discrete_covariance_ellipsoid([invrho(rho_disc*100), invrho(rho_disc*100)], D);
locmaxZ = simulateLocMax_Discrete_Stationary(D, sigma2, niter);
[empf, empz] = ecdf(locmaxZ);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z);

% generate another 1000 fields from same distribution for covariance estimation
for nn = 1:1000
    lat_data = normrnd( 0, 1, dim(1)+2*edge, dim(2)+2*edge);
    smoothed_fconv = fconv( lat_data, FWHM, 2);
    f_new(:,:,nn) = smoothed_fconv((edge+1):(edge+dim(1)),(edge+1):(edge+dim(2)));
end

f_new = f_new/std(f_new(:));
cov1 = empiricalCov(D,f_new);

% Make the covariance matrix spd
[V,DD] = eig(cov1);
DD = max(10e-10, DD);
sigma2 = V*DD*V';

locmaxZ3 = simulateLocMax_Discrete_Stationary(D, sigma2, niter);
[empf3, empz3] = ecdf(locmaxZ3);
empz3 = empz3(2:end); empf3 = empf3(2:end);   %remove non-unique point
pval_dist2 = 1-interp1(empz3,empf3,z);

figure();
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot(pval_ref,pval_dist,'LineWidth', 2);
xlim([0 0.05])
ylim([0 0.05])
hold on
plot(pval_ref,pval_dist2,'LineWidth', 2);
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
legend('Theoretical discrete covariance function', 'Empirical covariance function', '45 degree line', 'Location','northwest');

ax = gca;
set(gca,'FontSize', 18)
axis square
