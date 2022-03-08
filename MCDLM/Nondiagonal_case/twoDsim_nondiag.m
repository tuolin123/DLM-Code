%% % Description: This script provides an empirical distrbution of the 
% local maxima in a 2D gaussian random field by simulation from the 
% theoretical multivariate gaussian distribution.

%stddev_fwhm = 0.6;
%rho = round(exp(-1/2*1/(2*stddev_fwhm^2)),2);
D = 2;
sigma = 1;
niter = 1e6;
nondiag = 1;

%% Simulate the local maxima from the theoretical multivariate gaussian distribution
load('inverse_rho_discrete.mat')
rho_disc = 0.86;
stddev_fwhm = invrho(rho_disc*100);
rho = discrete_covariance(stddev_fwhm, D); % inverse of rho = 0.99, which should be the standard deviation 
tic
locmaxZ = simulateLocMax_Discrete(D,rho,sigma,niter,nondiag);
toc
%% %% Density of local maxima(continuous case)
kappa   = 1; % theoretical kappa
density = peakHeightDensity( 2, kappa );
dz = 0.001;
z=-2.5:dz:5; 
qcont = density(z);
%% %% Density of local maxima(discrete case)
nz=length(z);
ddrho = 0.001;
drho=0:ddrho:0.999;
nrho=length(drho);
p=zeros(nz,nrho);
nnb=ones(nz,3)*2;
q=ddlm(z,ones(nz,3)*rho_disc,nnb,2);
qnormalize=q./(sum(q)*dz);

%% %% Plot the density over simulation histogram
%histogram(locmaxZ, 'Normalization', 'pdf')
ng = 1000;
minz = min(locmaxZ);
maxz = max(locmaxZ);
pts = (minz-1):(maxz-minz+2)/ng:(maxz+1);
[g,xi] = ksdensity(locmaxZ, pts);
plot(xi, g)
hold on
plot(z,qcont,'g')
plot(z,qnormalize,'r')
legend({'simulation','continuous','discrete'},'Location', 'northwest')

%% %% get the localmax from the simulated field
locmaxZ2 = [];
nsim = 10000;
dim = [50 50];
f_new = zeros([dim(1:2),nsim]);
cut = 1;
FWHM = 2*sqrt(2*log(2))*stddev_fwhm;
edge = ceil(4*stddev_fwhm);

for nn = 1:nsim
    lat_data = normrnd( 0, 1, dim(1)+2*edge, dim(2)+2*edge);
    smoothed_fconv = fconv( lat_data, FWHM, 2);
    f_new(:,:,nn) = smoothed_fconv((edge+1):(edge+dim(1)),(edge+1):(edge+dim(2)));
end

f_new = f_new/std(f_new(:));

for nn = 1:nsim
    Z = f_new(:,:,nn);
    Z = squeeze(Z);
    %Z = interp3(Z, 'cubic');
    %Z = Z(range,range); % since two dimensional
    
    % Find local maxima of the field Z and remove maxima at the boundary
    % Imax = imregionalmax(Z); Imin = imregionalmin(Z);
    % Not including the diagonal
    Imax = imregionalmax(Z,4); Imin = imregionalmin(Z,4);
    
    % Not including the boundary
    Imax = Imax((1+cut):(end-cut), (1+cut):(end-cut));
    Imin = Imin((1+cut):(end-cut), (1+cut):(end-cut));
    
    % Rescale Z to match the dimensions of Imax and Imin
    Z = Z((1+cut):(end-cut), (1+cut):(end-cut));
    
    % add the new local maxima of Z to the vector containing all local
    % maxima
    locmaxZ2 = [locmaxZ2; Z(Imax); -Z(Imin)]; 
end

%% Get the smoothed eCDF of the simulated data
%pval_dist = pval_lookup(rho, z, D, "continuous", 0); % through the look-up table

[empf, empz] = ecdf(locmaxZ);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z);

% comparison of the look-up table method and simulated field method
figure();
%[empf, empz] = ecdf(locmaxZ2);
%empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
%pval_field = 1-interp1(empz,empf,z);
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot1 = plot(pval_ref,pval_dist,'LineWidth', 2);
plot1;
xlim([0 0.05])
ylim([0 0.05])
hold on
pval_cont = 1-cumsum(qcont).*dz;
plot2 = line(pval_ref,pval_cont, 'Color','r', 'LineWidth', 2);
plot2;
hold on
pval_disc = 1-cumsum(qnormalize).*dz;
plot3 = line(pval_ref,pval_disc, 'Color', 'green', 'LineWidth', 2);
plot3;
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
legend([plot2, plot3, plot1, hline], 'Theoretical continuous','Theoretical discrete','simulation (distribution)', '45 degree line', 'Location','northwest')
title(['$\rho$ = ' num2str(rho_disc)], 'Interpreter','latex')
ax = gca;
set(gca,'FontSize', 18)
axis square


