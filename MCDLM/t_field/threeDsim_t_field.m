%% % Description: This script provides an empirical distrbution of the 
% local maxima in a 2D gaussian random field by simulation from the 
% theoretical multivariate gaussian distribution.

load('inverse_rho_discrete.mat')
rho_disc = 0.5;
stddev_fwhm = invrho(rho_disc*100);
%rho = round(exp(-1/2*1/(2*stddev_fwhm^2)),2);
drho = 0.01;
rho = drho:drho:(1-drho);
D = 3;
sigma = 1;
niter = 1e5;
nu = 200;
nondiag = 0;

%% Simulate the local maxima from the theoretical multivariate gaussian distribution
tic
%locmaxZ = simulateLocMax(D,rho,sigma,niter,nondiag,nu);
rho_long = discrete_covariance(stddev_fwhm, D);
locmaxZ1 = [];
for j = 1:100
   locmaxZ = simulateLocMax_Discrete(D,rho_long,sigma,niter,nondiag,nu);
   locmaxZ1 = [locmaxZ locmaxZ1];
end
toc

%% %% Density of local maxima(continuous case)
kappa   = 1; % theoretical kappa
density = peakHeightDensity( D, kappa );
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

%% %% get the localmax from the simulated field
locmaxZ2 = [];
nsim = 100;
dim = [50 50 50];
f_new = zeros([dim(1:3),nsim]);
cut = 1;
FWHM = 2*sqrt(2*log(2))*stddev_fwhm;
edge = ceil(4*stddev_fwhm);

tic
for j = 1:10
    for nn = 1:nsim
        lat_data = normrnd( 0, 1, [dim(1)+2*edge, dim(2)+2*edge, dim(3)+2*edge, nu]);
        smoothed_fconv = convfield_t(lat_data, FWHM);
        f_new(:,:,:,nn) = smoothed_fconv.field((edge+1):(edge+dim(1)),(edge+1):(edge+dim(2)),(edge+1):(edge+dim(3)));
        
        Z = f_new(:,:,:,nn);
        Z = squeeze(Z);
        %Z = interp3(Z, 'cubic');
        %Z = Z(range,range); % since two dimensional

        % Find local maxima of the field Z and remove maxima at the boundary
        Imax = imregionalmax(Z); Imin = imregionalmin(Z);
        % Not including the diagonal
        % Imax = imregionalmax(Z,4); Imin = imregionalmin(Z,4);

        % Not including the boundary
        Imax = Imax((1+cut):(end-cut), (1+cut):(end-cut), (1+cut):(end-cut));
        Imin = Imin((1+cut):(end-cut), (1+cut):(end-cut), (1+cut):(end-cut));

        % Rescale Z to match the dimensions of Imax and Imin
        Z = Z((1+cut):(end-cut), (1+cut):(end-cut), (1+cut):(end-cut));

        % add the new local maxima of Z to the vector containing all local
        % maxima
        locmaxZ2 = [locmaxZ2; Z(Imax); -Z(Imin)]; 
    end
end
toc

%% Get the smoothed eCDF of the simulated data
%pval_dist = pval_lookup(rho, z, D, "continuous", 1); % through the look-up table

[empf, empz] = ecdf(locmaxZ1);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z); % without lookup table

% comparison of the look-up table method and simulated field method
figure();
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
set(gca,'FontSize', 18)
xlabel('$p$-value', 'Interpreter', 'latex', 'fontsize',24)
ylabel('ecdf($p$)', 'Interpreter', 'latex', 'fontsize',24)
legend([plot2, plot1, plot3, hline], 'Theoretical continuous','Fully connected DLM','Partially connected DLM', '45 degree line', 'Location','northwest')
title(['$\rho$ = ' num2str(rho_disc) '\ (' 'FWHM = ' num2str(round(FWHM,1)) ')'], 'Interpreter','latex', 'fontsize',24)
%ax = gca;
axis square

