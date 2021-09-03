mask = spm_read_vols(spm_vol('MNImask.nii'));

mask_2d = squeeze(mask(:,:,45));
mask_2d_erode = dilate_mask(mask_2d, -1);

Dim = [91,109,91];
nimgs = 198;

%% Use grey and white mask only (separate)
fid1 = fopen('avg152T1_white.img'); 
fid2 = fopen('avg152T1_gray.img');
data = fread(fid1);
GM_intensity = reshape(data, [91,109,91])/255;
%imagesc(GM_intensity(:,:,45));
colorbar

mask = GM_intensity > 0.5;
imagesc(mask(:,:,45))

% use the original images with the white matter mask
mask_2d = squeeze(mask(:,:,45));
mask_2d_erode = dilate_mask(mask_2d, -1);

%% try to use white noise with mask
noise = zeros([Dim(1:2),nimgs]);
mu = 0; sigma = 1;
for i = 1:nimgs
    noise(:,:,i) = mu+sigma*randn(Dim(1:2));
end

imgs_white = zeros([Dim(1:2) nimgs]);
for I = 1:nimgs
    imgs_white(:,:,I) = squeeze(noise(:,:,I)).*mask_2d;
    %imgs_white(:,:,I) = squeeze(noise(:,:,I));
end

%% Try Scale 3 Laplacian field 
whitefield = wfield(Dim(1:2), nimgs, 'L', 3);
noise = whitefield.field;
imgs_white = zeros([Dim(1:2) nimgs]);
for I = 1:nimgs
    imgs_white(:,:,I) = squeeze(noise(:,:,I)).*mask_2d;
end

imgs_white = zero2nan(imgs_white);
f_new = imgs_white/std(imgs_white(~isnan(imgs_white)));

%% Gaussian noise smoothed with a kernel
f_new = zeros([Dim(1:2), nimgs]);
cut = 1;
stddev_fwhm = 0.6;
FWHM = 2*sqrt(2*log(2))*stddev_fwhm;
edge = ceil(4*stddev_fwhm); 

for nn = 1:nimgs
    lat_data = normrnd( 0, 1, Dim(1)+2*edge, Dim(2)+2*edge);
    smoothed_fconv = fconv( lat_data, FWHM, 2);
    f_new(:,:,nn) = smoothed_fconv((edge+1):(edge+Dim(1)),(edge+1):(edge+Dim(2)));
    f_new(:,:,nn) = squeeze(f_new(:,:,nn)).*zero2nan(mask_2d);
end

%f_new = f_new/std(f_new(~isnan(f_new))); % standardize voxelwise
f_new = f_new./std(f_new,0,3);

%% Try the real data (use imgs_cat from "beijing.m")
f_new = sqrt(max(0,squeeze(imgs_cat(:,:,45,:)))) - sqrt(max(0,-squeeze(imgs_cat(:,:,45,:))));
%f_new = squeeze(imgs_cat(:,:,45,:));
f_new = (f_new-mean(f_new, D+1)).*zero2nan(mask_2d);
f_new = f_new./std(f_new,0,3);
%f_new = -f_new./std(f_new,0,3);


%% Get the localmax
locmaxZ2 = [];
locmaxZ = [];
B = 1000;
nb = 40;
sigma = 1;
imgs_new = ones([Dim(1) Dim(2) B]);
sigma_all = zeros(1,B);

tic
for k = 1:B
    samp_ind = randsample(1:nimgs, nb);
    imgs_mask_2d_samp = f_new(:,:,samp_ind);
    Z = sqrt(nb)*mean(imgs_mask_2d_samp,3); %take the average and multiply by root-n
    %Z = imgs_white(:,:,k);
    [peaklocs, peakinds, peakvals] = lmindices(Z,'all',mask_2d);
    Z2 = zeros(Dim(1:2));
    Z2(peakinds) = peakvals;
    Z2 = Z2.*mask_2d_erode;
    peakvals = Z2(Z2 ~= 0);
    locmaxZ2 = [locmaxZ2;peakvals];
    
    % estimate the rho % 
    [rho_est, sigma_est] = estsd(Z, mask_2d, D, 0.5);
    sigma_all(k) = mean(sigma_est);
    
    % sigma1 = discrete_covariance_ellipsoid(sigma_est, 2);
    rho = discrete_covariance(mean(sigma_est),D);
    niter = 1e5;
    %n_lim = 1e5;
    n_lim = 1e3;

    locmaxZ1 = [];
    nj = 10000;
    
    for j = 1:nj
        %locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(D-1, sigma1, niter);
        locmaxZ_sigma1 = simulateLocMax_Discrete(D,rho,sigma,niter);
        locmaxZ1 = [locmaxZ_sigma1 locmaxZ1];
        if length(locmaxZ1) > n_lim
            break;
        end
    end

    locmaxZ = [locmaxZ1 locmaxZ];
    
    % create the new sample
    imgs_new(:,:,k) = Z; 
end
toc

%% %% Density of local maxima(continuous case)
kappa   = 1; % theoretical kappa
density = peakHeightDensity( length(Dim)-1, kappa );
dz = 0.001;
z = min(locmaxZ):dz:max(locmaxZ); 
qcont = density(z);
%% %% Density of local maxima(discrete case)
sigma_est = estsd(f_new, mask_2d, D);
rho_disc = exp(-1/2*1./(2*sigma_est.^2));

nz=length(z);
ddrho = 0.001;
drho=0:ddrho:0.999;
nrho=length(drho);
p=zeros(nz,nrho);
nnb=ones(nz,D-1)*2;
q=ddlm(z,ones(nz,D-1).*rho_disc,nnb,D-1);
qnormalize=q./(sum(q)*dz);

%% plot
[empf, empz] = ecdf(locmaxZ);
empz = empz(2:end); empf = empf(2:end);   %remove non-unique point
pval_dist = 1-interp1(empz,empf,z);

figure();
[ksf, ksz] = ksdensity(locmaxZ2, z);
pval_ref = 1-cumsum(ksf)*dz;

plot(pval_dist,pval_ref,'LineWidth', 2);
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
hold on
hline = refline(1,0);
set(hline,'LineStyle',':');
set(hline,'Color','black');
set(hline,'LineWidth', 2);
set(gca,'FontSize', 18)
xlabel('$p$-value', 'Interpreter', 'latex', 'fontsize',24)
ylabel('Reference $p$-value', 'Interpreter', 'latex', 'fontsize',24)
legend('Fully connected DLM','Continuous','Partially connected DLM',...
     '45 degree line', 'Location','northwest')
title("Beijing 4mm")
ax = gca;

axis square

%% check histogram
% after square-root
histogram(imgs_new(~isnan(imgs_new)), 'Normalization', 'PDF')
hold on
x = -5:0.01:5;
y = normpdf(x);
plot(x,y)
legend('sqrt rsd','N(0,1)','Location','northwest')

% before square-root
histogram(imgs_new(~isnan(imgs_new)), 'Normalization', 'PDF')
hold on
x = -5:0.01:5;
y = normpdf(x);
plot(x,y)
legend('rsd','N(0,1)','Location','northwest')

qqplot(imgs_new(~isnan(imgs_new)))

