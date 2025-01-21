%% Average of two ellipsoid example
% Generate the field
load('inverse_rho_discrete.mat')

rho_disc1 = 0.01;
rho_disc2 = 0.5;

stddev_fwhm1 = rho2sd(rho_disc1);
stddev_fwhm2 = rho2sd(rho_disc2);
Dim = 2;

nsim = 100;
dim = [50 50];
f_new = zeros([dim(1:2),nsim]);
f_new1 = zeros([dim(1:2),nsim]);
f_ave = zeros([dim(1:2),nsim]);
cut = 1;
FWHM_1 = 2*sqrt(2*log(2))*stddev_fwhm1;
FWHM_2 = 2*sqrt(2*log(2))*stddev_fwhm2;
edge1 = ceil(4*stddev_fwhm1);
edge2 = ceil(4*stddev_fwhm2);

% for nn = 1:nsim
%     lat_data = wfield([dim(1)+2*edge1, dim(2)+2*edge2]); resadd = 0; enlarge = 0;
%     params = ConvFieldParams([FWHM_1, FWHM_2], resadd, enlarge);
%     smoothed_fconv = convfield(lat_data, params);
%     f_new(:,:,nn) = smoothed_fconv.field((edge1+1):(edge1+dim(1)),(edge2+1):(edge2+dim(2)));
% end
% 
% f_new = f_new/std(f_new(:));
% %imagesc(f_new(:,:,1))
% 
% for nn = 1:nsim
%     lat_data = wfield([dim(2)+2*edge2, dim(1)+2*edge1]); resadd = 0; enlarge = 0;
%     params = ConvFieldParams([FWHM_2, FWHM_1], resadd, enlarge);
%     smoothed_fconv = convfield(lat_data, params);
%     f_new1(:,:,nn) = smoothed_fconv.field((edge2+1):(edge2+dim(2)),(edge1+1):(edge1+dim(1)));
% end
% 
% f_new1 = f_new1/std(f_new1(:));
% 
% for nn = 1:nsim
%     f_ave(:,:,nn) = (f_new(:,:,nn) + f_new1(:,:,nn))/2;
% end

%% Simulate locmax from the gaussian distribution
sigma = discrete_covariance_ellipsoid([stddev_fwhm1 stddev_fwhm2], 2);
sigma_flip = discrete_covariance_ellipsoid([stddev_fwhm2 stddev_fwhm1], 2);
sigma1 = (sigma + sigma_flip)/4; 
niter = 1e5; 
n_lim = 1e6;

locmaxZ1 = [];
nj = 10000;

for j = 1:nj
    locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(Dim, sigma1, niter);
    locmaxZ1 = [locmaxZ_sigma1 locmaxZ1];
    if length(locmaxZ1) > n_lim
        break;
    end
end

% Estimatating the empirical covariance method
% locmaxZ = [];
% tic
% sigma2 = empiricalCov(Dim,f_new);
% 
% % Make the covariance matrix spd
% [V,D] = eig(sigma2);
% D = max(10e-10, D);
% sigma2 = V*D*V';
% for j = 1:nj
%     locmaxZ_sigma2 = simulateLocMax_Discrete_Stationary(Dim, sigma2, niter);
%     locmaxZ = [locmaxZ_sigma2 locmaxZ];
%     if length(locmaxZ) > n_lim
%         break;
%     end
% end
% toc

%% Simulate locmax from the field
tic
locmaxZ2 = [];
for j = 1:50
    for nn = 1:nsim
        lat_data = wfield([dim(1)+2*edge1, dim(2)+2*edge2]); resadd = 0; enlarge = 0;
        params = ConvFieldParams([FWHM_1, FWHM_2], resadd, enlarge);
        smoothed_fconv = convfield(lat_data, params);
        f_new(:,:,nn) = smoothed_fconv.field((edge1+1):(edge1+dim(1)),(edge2+1):(edge2+dim(2)));
    end

    f_new = f_new/std(f_new(:));

    
    for nn = 1:nsim
        lat_data = wfield([dim(2)+2*edge2, dim(1)+2*edge1]); resadd = 0; enlarge = 0;
        params = ConvFieldParams([FWHM_2, FWHM_1], resadd, enlarge);
        smoothed_fconv = convfield(lat_data, params);
        f_new1(:,:,nn) = smoothed_fconv.field((edge2+1):(edge2+dim(2)),(edge1+1):(edge1+dim(1)));
    end

    f_new1 = f_new1/std(f_new1(:));

    for nn = 1:nsim
        Z = (f_new(:,:,nn) + f_new1(:,:,nn))/2;
        Z = squeeze(Z);
        %Z = interp3(Z, 'cubic');
        %Z = Z(range,range); % since two dimensional

        % Find local maxima of the field Z and remove maxima at the boundary
        Imax = imregionalmax(Z); Imin = imregionalmin(Z);
        % Not including the diagonal
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
toc

%% Histogram comparison

histogram(locmaxZ2, 'Normalization', 'pdf')
hold on
histogram(locmaxZ1, 'Normalization', 'pdf')

%% Figure
dz = 0.001;
z=-2.5:dz:5; 

[empf1, empz1] = ecdf(locmaxZ1);
empz1 = empz1(2:end); empf1 = empf1(2:end);   %remove non-unique point
pval_dist1 = 1-interp1(empz1,empf1,z);

figure();
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot(pval_dist1,pval_ref,'LineWidth', 2);
xlim([0 0.05])
ylim([0 0.05])
hold on
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
set(gca,'FontSize', 18)
xlabel('$p$-value', 'Interpreter', 'latex', 'fontsize',22)
ylabel('Reference $p$-value', 'Interpreter', 'latex', 'fontsize',22)
yticks(0:0.01:0.05)
lgd = legend('MCDLM', '45 degree line', 'Location','southeast');
lgd.FontSize = 18;
title({['$\rho_1$ = ' num2str(rho_disc1) '\ (' 'FWHM = ' num2str(round(FWHM_1,1)) ')']...
    ['$\rho_2$ = ' num2str(rho_disc2) '\ (' 'FWHM = ' num2str(round(FWHM_2,1)) ')']},...
    'Interpreter','latex','fontsize',22)
ax = gca;

axis square

saveas(gcf, '2D_nonsep_rho1_0.01_rho2_0.5.jpg')