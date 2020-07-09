%% % Description: This script provides an empirical distrbution of the 
% local maxima in a 2D gaussian random field by simulation from the 
% theoretical multivariate gaussian distribution.

stddev = 5.1;
rho = round(exp(-1/2*1/(2*stddev^2)),2);
%rho = exp(-1/2*1/(2*stddev^2));
D = 2;
sigma = 1;
niter = 1e5;
nondiag = 1;

%% Simulate the local maxima from the theoretical multivariate gaussian distribution
tic
locmaxZ = simulateLocMax(D,rho,sigma,niter,nondiag);
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
q=ddlm(z,ones(nz,3)*rho,nnb,2);
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

%% %% Plot the histogram of local maxima from simulated field
locmaxZ2 = [];
path_data = 'isoL505030nsim10000n1_gauss_stddev5.1.mat';
load( path_data );
slice = 1;

% Loop over realisations to find the maxima in each field
tic
for nn = 1:nsim
    % Get the sample field
    Z = f(:,:,slice,nn);
    Z = squeeze(Z);
    %Z = interp3(Z, 'cubic');
    
    % Find local maxima of the field Z and remove maxima at the boundary
    % Not including the diagonal
    Imax = imregionalmax(Z,4); Imin = imregionalmin(Z,4);
    Imax = Imax( (1+cut):(end-cut), (1+cut):(end-cut));
    Imin = Imin((1+cut):(end-cut), (1+cut):(end-cut));
    
    % Rescale Z to match the dimensions of Imax and Imin
    Z = Z((1+cut):(end-cut), (1+cut):(end-cut));
    
    % add the new local maxima of Z to the vector containing all local
    % maxima
    locmaxZ2 = [locmaxZ2; Z(Imax); -Z(Imin)]; 
    %locmaxZ2 = [locmaxZ2; Z(Imax)];
end
toc
%histogram(locmaxZ2, 'Normalization', 'pdf')
%hold on

%% Get the smoothed eCDF of the simulated data
load('presmooth_lookup_2D_nondiag.mat')
load('lookup_2D_CDF_nondiag.mat')
dp = cdf_table2D.p;
dx = cdf_table2D.x;
drho = cdf_table2D.rho;
dp(dp<0) = 0;
[X,Y] = meshgrid(dx,drho);

pval_dist = 1-interp2(X,Y,dp,z,rho);
pval_dist_pre = 1-interp2(X,Y,matcdf,z,rho);

% comparison of the look-up table method and simulated field method
figure();
%[empf, empz] = ecdf(locmaxZ2);
%empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
%pval_field = 1-interp1(empz,empf,z);
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot(pval_ref,pval_dist,'LineWidth', 2);
xlim([0 0.05])
ylim([0 0.05])
hold on
pval_cont = 1-cumsum(qcont).*dz;
line(pval_ref,pval_cont, 'Color','r', 'LineWidth', 2);
hold on
pval_disc = 1-cumsum(qnormalize).*dz;
line(pval_ref,pval_disc, 'Color', 'green', 'LineWidth', 2);
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
legend('Look-up table','Theoretical continuous', 'Theoretical discrete', '45 degree line', 'Location','northwest')
title(['$\rho$ = ' num2str(rho)], 'Interpreter','latex')
ax = gca;
set(gca,'FontSize', 18)
axis square


