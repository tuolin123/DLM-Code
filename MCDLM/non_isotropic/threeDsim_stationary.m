%% Ellipsoid example
% Generate the field
load('inverse_rho_discrete.mat')

rho_disc1 = 0.95;
rho_disc2 = 0.95;
rho_disc3 = 0.95;

stddev_fwhm1 = invrho(rho_disc1*100);
stddev_fwhm2 = invrho(rho_disc2*100);
stddev_fwhm3 = invrho(rho_disc3*100);
Dim = 3;

nsim = 200;
dim = [50 50 50];
f_new = zeros([dim(1:3),nsim]);
cut = 1;
FWHM_1 = 2*sqrt(2*log(2))*stddev_fwhm1;
FWHM_2 = 2*sqrt(2*log(2))*stddev_fwhm2;
FWHM_3 = 2*sqrt(2*log(2))*stddev_fwhm3;
edge1 = ceil(4*stddev_fwhm1);
edge2 = ceil(4*stddev_fwhm2);
edge3 = ceil(4*stddev_fwhm3);

for nn = 1:nsim
    lat_data = wfield([dim(1)+2*edge1, dim(2)+2*edge2, dim(3)+2*edge3]); resadd = 0; enlarge = 0;
    params = ConvFieldParams([FWHM_1, FWHM_2, FWHM_3], resadd, enlarge);
    smoothed_fconv = convfield(lat_data, params);
    f_new(:,:,:,nn) = smoothed_fconv.field((edge1+1):(edge1+dim(1)),(edge2+1):(edge2+dim(2)), (edge3+1):(edge3+dim(2)));
end

f_new = f_new/std(f_new(:));
%imagesc(f_new(:,:,1))

%% Simulate locmax from the gaussian distribution
sigma1 = discrete_covariance_ellipsoid([stddev_fwhm1 stddev_fwhm2, stddev_fwhm3], Dim);
% [V,D] = eig(sigma1);
% D = max(10e-10, D);
% sigma1 = V*D*V';
niter = 1e5;

if min([rho_disc1, rho_disc2, rho_disc3]) == 0.01
    n_lim = 1e6;
elseif min([rho_disc1, rho_disc2, rho_disc3]) == 0.5
    n_lim = 1e6;
elseif min([rho_disc1, rho_disc2, rho_disc3]) == 0.99
    n_lim = 2e5;
end

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
locmaxZ = [];
tic
sigma2 = empiricalCov(Dim,f_new);

% Make the covariance matrix spd
[V,D] = eig(sigma2);
D = max(10e-10, D);
sigma2 = V*D*V';
for j = 1:nj
    locmaxZ_sigma2 = simulateLocMax_Discrete_Stationary(Dim, sigma2, niter);
    locmaxZ = [locmaxZ_sigma2 locmaxZ];
    if length(locmaxZ) > n_lim
        break;
    end
end
toc

%% Simulate locmax from the field
locmaxZ2 = [];
nsim = 1000;
for j = 1:5
    for nn = 1:nsim
        lat_data = wfield([dim(1)+2*edge1, dim(2)+2*edge2, dim(3)+2*edge3]); resadd = 0; enlarge = 0;
        params = ConvFieldParams([FWHM_1, FWHM_2, FWHM_3], resadd, enlarge);
        smoothed_fconv = convfield(lat_data, params);
        f_new(:,:,:,nn) = smoothed_fconv.field((edge1+1):(edge1+dim(1)),(edge2+1):(edge2+dim(2)),(edge3+1):(edge3+dim(3)));
    end

    f_new = f_new/std(f_new(:));

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

%% Histogram comparison

histogram(locmaxZ2, 'Normalization', 'pdf')
hold on
histogram(locmaxZ, 'Normalization', 'pdf')

%% combining plots together
load(strcat('3D_rho1_',num2str(rho_disc1),'_rho2_',num2str(rho_disc2),'_nsim_200.mat'));
% locmaxZ3 = load(strcat('3D_rho1_',num2str(rho_disc1),'_rho2_',num2str(rho_disc2),'_nsim_1000.mat'));
% locmaxZ3 = locmaxZ3.locmaxZ;
%locmaxZ4 = load(strcat('3D_rho1_',num2str(rho_disc1),'_rho2_',num2str(rho_disc2),'_nsim_10000.mat'));
%locmaxZ4 = locmaxZ4.locmaxZ;
locmaxZ5 = load(strcat('3D_rho1_',num2str(rho_disc1),'_rho2_',num2str(rho_disc2),'_nsim_20.mat'));
locmaxZ5 = locmaxZ5.locmaxZ;
locmaxZ6 = load(strcat('3D_rho1_',num2str(rho_disc1),'_rho2_',num2str(rho_disc2),'_nsim_50.mat'));
locmaxZ6 = locmaxZ6.locmaxZ;
locmaxZ7 = load(strcat('3D_rho1_',num2str(rho_disc1),'_rho2_',num2str(rho_disc2),'_nsim_100.mat'));
locmaxZ7 = locmaxZ7.locmaxZ;

dz = 0.001;
z=-2.5:dz:5; 

[empf, empz] = ecdf(locmaxZ);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z);

[empf1, empz1] = ecdf(locmaxZ1);
empz1 = empz1(2:end); empf1 = empf1(2:end);   %remove non-unique point
pval_dist1 = 1-interp1(empz1,empf1,z);

% [empf2, empz2] = ecdf(locmaxZ3);
% empz2 = empz2(2:end); empf2 = empf2(2:end);   %remove non-unique point
% pval_dist2 = 1-interp1(empz2,empf2,z);

% [empf3, empz3] = ecdf(locmaxZ4);
% empz3 = empz3(2:end); empf3 = empf3(2:end);   %remove non-unique point
% pval_dist3 = 1-interp1(empz3,empf3,z);

[empf4, empz4] = ecdf(locmaxZ5);
empz4 = empz4(2:end); empf4 = empf4(2:end);   %remove non-unique point
pval_dist4 = 1-interp1(empz4,empf4,z);

[empf5, empz5] = ecdf(locmaxZ6);
empz5 = empz5(2:end); empf5 = empf5(2:end);   %remove non-unique point
pval_dist5 = 1-interp1(empz5,empf5,z);

[empf6, empz6] = ecdf(locmaxZ7);
empz6 = empz6(2:end); empf6 = empf6(2:end);   %remove non-unique point
pval_dist6 = 1-interp1(empz6,empf6,z);

figure();
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot(pval_dist1,pval_ref,'LineWidth', 2);
xlim([0 0.05])
ylim([0 0.05])
hold on
plot(pval_dist4,pval_ref,'LineWidth', 2);
hold on
plot(pval_dist5,pval_ref,'LineWidth', 2);
hold on
plot(pval_dist6,pval_ref,'LineWidth', 2);
hold on
plot(pval_dist,pval_ref,'LineWidth', 2);
hold on
% plot(pval_dist2,pval_ref,'LineWidth', 2);
% hold on
%plot(pval_dist3,pval_ref,'LineWidth', 2);
%hold on
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
set(gca,'FontSize', 15)
xlabel('$p$-value', 'Interpreter', 'latex', 'fontsize',18)
ylabel('Reference $p$-value', 'Interpreter', 'latex', 'fontsize',18)
xticks(0:0.01:0.05)
lgd = legend('Tcf', 'Ecf (nsim = 20)','Ecf (nsim = 50)',...
    'Ecf (nsim = 100)','Ecf (nsim = 200)',...
    '45 degree line', 'Location','southeast');
lgd.FontSize = 10;
title({['$\rho_1$ = $\rho_2$ = $\rho_3$ = ' num2str(rho_disc1) '\ (' 'FWHM = ' num2str(round(FWHM_1,1)) ')']},...
    'Interpreter','latex')
ax = gca;

axis square
%saveas(gcf, '3D_rho1_0.01_rho2_0.01_rho3_0.01.jpg')