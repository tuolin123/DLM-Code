%% % Description: This script provides an empirical distrbution of the 
% local maxima in a 3D gaussian random field by simulation from the 
% theoretical multivariate gaussian distribution.

D = 3;
sigma = 1;
niter = 1e5;

%% Simulate the local maxima from the theoretical multivariate gaussian distribution
load('inverse_rho_discrete.mat')
rho_disc = 0.99;
stddev_fwhm = invrho(rho_disc*100);
if rho_disc == 0.01
    n_lim = 1e6;
elseif rho_disc == 0.5
    n_lim = 1e6;
elseif rho_disc == 0.99
    n_lim = 2e5;
end

tic
locmaxZ1 = [];
rho = discrete_covariance(stddev_fwhm, D); % inverse of rho = 0.99, which should be the standard deviation 
for j = 1:10000
    locmaxZ = simulateLocMax_Discrete(D,rho,sigma,niter);
    locmaxZ1 = [locmaxZ locmaxZ1];
    if length(locmaxZ1) > n_lim
        break;
    end
end
toc

%% After this find rho corresponding to specific standard deviation
%stddev_fwhm = 5;
%[M, I] = min(abs(invrho-stddev_fwhm));
%rho_disc = rho(I);
%rho_disc = 0.01;

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
q=ddlm(z,ones(nz,3)*rho_disc,nnb,3);
qnormalize=q./(sum(q)*dz);

%% %% get the localmax from the simulated field
locmaxZ2 = [];
nsim = 5000;
dim = [50 50 50];
f_new = zeros([dim,nsim]);
cut = 1;
FWHM = 2*sqrt(2*log(2))*stddev_fwhm;
edge = ceil(4*stddev_fwhm);

tic
for i = 1:2
    for nn = 1:nsim
        lat_data = normrnd( 0, 1, dim(1)+2*edge, dim(2)+2*edge, dim(3)+2*edge);
        smoothed_fconv = fconv( lat_data, FWHM, 3);
        f_new(:,:,:,nn) = smoothed_fconv((edge+1):(edge+dim(1)),(edge+1):(edge+dim(2)),(edge+1):(edge+dim(3)));
    end

    f_new = f_new/std(f_new(:));
    % save(['iso',num2str(dim(1)),num2str(dim(2)),num2str(dim(3)),'nsim'...
    %                 ,num2str(nsim),'_stddev',num2str(stddev_fwhm),'.mat'])

    for nn = 1:nsim
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
tic
pval_dist = pval_lookup(rho_disc, z, D, 'discrete',1); % through the look-up table
toc
% [empf, empz] = ecdf(locmaxZ1);
% empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
% pval_dist = 1-interp1(empz,empf,z);

% comparison of the look-up table method and simulated field method
figure();
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

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
lgd = legend([plot3, plot1, plot2, hline], 'ADLM', 'MCDLM', 'Continuous RFT','45 degree line', 'Location','southeast');
lgd.FontSize = 12;
%legend([plot2, plot1, plot3, hline], 'Theoretical continuous','Look-up table','Partially connected DLM', '45 degree line', 'Location','southeast')
title(['$\rho$ = ' num2str(rho_disc) '\ (' 'FWHM = ' num2str(round(FWHM,1)) ')'], 'Interpreter','latex', 'fontsize',18)
%ax = gca;
axis square

saveas(gcf, 'lookup_3D_rho0.99_disc.jpg')
%% Check the independent case
rhoindp = discrete_covariance(0, D);
locmaxZ = simulateLocMax_Discrete(D,rhoindp,sigma,niter);
[find, zind] = ecdf(locmaxZ);
zind = zind(2:end); find = find(2:end);   %remove non-unique point
pval_dist = 1-interp1(zind,find,z);

for nn = 1:nsim
    lat_data = normrnd( 0, 1, dim(1), dim(2));
    f_new(:,:,nn) = lat_data;
    
    Z = f_new(:,:,nn);
    Z = squeeze(Z);
    %Z = interp3(Z, 'cubic');
    %Z = Z(range,range); % since two dimensional
    
    % Find local maxima of the field Z and remove maxima at the boundary
    Imax = imregionalmax(Z); Imin = imregionalmin(Z);
    % Not including the diagonal
    % Imax = imregionalmax(Z,8); Imin = imregionalmin(Z,8);
    Imax = Imax( (1+cut):(end-cut), (1+cut):(end-cut));
    Imin = Imin((1+cut):(end-cut), (1+cut):(end-cut));
    
    % Rescale Z to match the dimensions of Imax and Imin
    Z = Z((1+cut):(end-cut), (1+cut):(end-cut));
    
    % add the new local maxima of Z to the vector containing all local
    % maxima
    locmaxZ2 = [locmaxZ2; Z(Imax); -Z(Imin)]; 
end
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;
plot(pval_ref,pval_dist,'LineWidth', 2);
xlim([0 0.05])
ylim([0 0.05])
hold on
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
legend('Discrete covaraince function', '45 degree line', 'Location','northwest')
title(['$\rho$ = ' num2str(0)], 'Interpreter','latex')
ax = gca;
set(gca,'FontSize', 18)
axis square

%% Comparing with continuous covariance function
rho_vec = discrete_covariance(invrho(ceil(100*rho_disc)), D);
locmaxZ = simulateLocMax_Discrete(D,rho_vec,sigma,niter);
[empf, empz] = ecdf(locmaxZ);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z); % without lookup table

rho_cont = round(exp(-1/2*1/(2*stddev_fwhm^2)),2); % continuous rho
locmaxZ3 = simulateLocMax(D,rho_cont,sigma,niter); % continuous covariance function
[empf3, empz3] = ecdf(locmaxZ3);
empz3 = empz3(2:end); empf3 = empf3(2:end);   %remove non-unique point
pval_dist2 = 1-interp1(empz3,empf3,z);

nondiag = 1;
locmaxZ4 = simulateLocMax(D,rho_cont,sigma,niter,nondiag); % continuous covariance function
[empf4, empz4] = ecdf(locmaxZ4);
empz4 = empz4(2:end); empf4 = empf4(2:end);   %remove non-unique point
pval_dist3 = 1-interp1(empz4,empf4,z);

% comparison of the look-up table method and simulated field method
figure();
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot1 = plot(pval_ref,pval_dist,'LineWidth', 2);
xlim([0 0.05])
ylim([0 0.05])
hold on
plot2 = plot(pval_ref,pval_dist2,'LineWidth', 2);
hold on
plot3 = plot(pval_ref,pval_dist3,'LineWidth', 2);
hold on
pval_cont = 1-cumsum(qcont).*dz;
plot4 = line(pval_ref,pval_cont, 'Color','r', 'LineWidth', 2);
hold on
pval_disc = 1-cumsum(qnormalize).*dz;
plot5 = line(pval_ref,pval_disc, 'Color', 'green', 'LineWidth', 2);
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
[hleg, hobj] = legend([plot4, plot2, plot1, plot3, plot5, hline],...
        'Theoretical continuous','Continuous covariance function',...
        'Discrete covaraince function','continuous covariance function nondiag',...
        'Theoretical discrete', '45 degree line', 'Location','northwest');
set(hleg,'Location','northwest','FontSize',7);
title(['$\rho$ = ' num2str(rho_disc)], 'Interpreter','latex')
ax = gca;
set(gca,'FontSize', 18)
axis square

