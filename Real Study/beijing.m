%% read in data
% try the 4mm smoothing kernel one
imagefiles = dir('Beijing*boxcar10*4mm*gz'); 
nfiles = length(imagefiles);
imgs = cell(nfiles);

for ii=1:nfiles
   filename = imagefiles(ii).name;
   img = spm_read_vols(spm_vol(filename));
   imgs{ii} = img;
end

nimgs = length(imgs);
Dim = size(imgs{1});
D = length(Dim);

% concatenate all the cells together into one array
imgs_cat = zeros([Dim nimgs]);
for i = 1:nimgs
    imgs_cat(:,:,:,i) = imgs{i};
end

%% standardize and apply the mask
mask = spm_read_vols(spm_vol('MNImask.nii'));
mask_na = zero2nan(mask);

%masked_image = cellfun(@(x) x.*mask, imgs, 'UniformOutput', false);

%edit lmindicies %% top = "all"

% check the image
% figure;
% subplot(2,1,1)
% imagesc(imgs{1}(:,:,50))
% 
% subplot(2,1,2)
% imagesc(masked_image{1}(:,:,50))

% masked_img_cat = zeros([Dim nimgs]);
% for i = 1:nimgs
%     masked_img_cat(:,:,:,i) = masked_image{i}/std(masked_image{i}(:));
% end

% mean
mean_over_subjects = mean(imgs_cat, length(Dim)+1);

imgs_mask = zeros([Dim nimgs]);

% Subtract the mean and multiply by the mask.
for I = 1:nimgs
    imgs_mask(:,:,:,I) = (imgs_cat(:,:,:,I) - mean_over_subjects).*mask_na;
end

% standardize for each voxel across all images
std_dev = std(imgs_mask, 0, 4);
imgs_mask = imgs_mask./std_dev;

% check the std
std(imgs_mask(~isnan(imgs_mask)))

%% estimate the FWHM (and nu for smoothing kernel), need RFTtoolbox
% nsubj = 20;
% 
% samp_ind = randsample(1:length(imgs), nsubj);
% imgs_samp = imgs(samp_ind);
% imgs_samp_cat = zeros([Dim nsubj]);
% for i = 1:nsubj
%     imgs_samp_cat(:,:,:,i) = imgs_samp{i};
% end
% 
% % check the figure to make sure
% figure;
% subplot(2,1,1)
% imagesc(imgs_samp{1}(:,:,50))
% 
% subplot(2,1,2)
% imagesc(imgs_samp_cat(:,:,50,1))
% 
% [fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est] = est_smooth(imgs_samp_cat, mask);
% 
% masked_imgsamp = masked_image(samp_ind); 
% masked_imgsamp_cat = zeros([Dim nsubj]);
% for i = 1:nsubj
%     masked_imgsamp_cat(:,:,:,i) = masked_imgsamp{i}/std(masked_imgsamp{i}(:));
% end
% covv = empiricalCov(3,masked_imgsamp_cat);
% covv1 = covv/covv(1,1);

%% Get the localmax
locmaxZ2 = [];
locmaxZ = [];
imgs_mask_2d = squeeze(imgs_mask(:,:,45,:));
mask_2d = squeeze(mask(:,:,45));
B = 200;
nb = 40;
imgs_new = ones([Dim(1) Dim(2) B]);
load('inverse_rho_discrete.mat')

tic
for k = 1:B
    samp_ind = randsample(1:nimgs, nb);
    imgs_mask_2d_samp = imgs_mask_2d(:,:,samp_ind);
    Z = sqrt(nb)*mean(imgs_mask_2d_samp,3); %take the average and multiply by root-n
    [peaklocs, peakinds, peakvals] = lmindices(Z,"all",mask_2d);
    locmaxZ2 = [locmaxZ2;peakvals];
    
    % estimate the rho
    Z_stand = Z/std(Z(~isnan(Z)));
    Z_mean = mean(Z(~isnan(Z)));
    rho_hor = (Z_stand(1:(Dim(1)-1),1:Dim(2))-Z_mean).*(Z_stand(2:Dim(1),1:Dim(2))-Z_mean);
    rho_vert = (Z_stand(1:Dim(1),1:(Dim(2)-1))-Z_mean).*(Z_stand(1:Dim(1),2:Dim(2))-Z_mean);
    rho = [mean(rho_hor(~isnan(rho_hor))) mean(rho_vert(~isnan(rho_vert)))];
    sigma_est = [invrho(round(rho(1),2)*100) invrho(round(rho(2),2)*100)];
    
    sigma1 = discrete_covariance_ellipsoid(sigma_est, D-1);
    niter = 1e5;
    n_lim = 1e5;

    locmaxZ1 = [];
    nj = 10000;

    for j = 1:nj
        locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(D-1, sigma1, niter);
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


%% Use full-connectivity DLM method
locmaxZ = [];
tic

[fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est] = est_smooth(imgs_new, mask_2d);

sigma1 = discrete_covariance_ellipsoid(sigma_est, D-1);
niter = 1e5;
n_lim = 1e6;

locmaxZ1 = [];
nj = 10000;

for j = 1:nj
    locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(D-1, sigma1, niter);
    locmaxZ1 = [locmaxZ_sigma1 locmaxZ1];
    if length(locmaxZ1) > n_lim
        break;
    end
end

locmaxZ = [locmaxZ1 locmaxZ];

toc

%% %% Density of local maxima(continuous case)
kappa   = 1; % theoretical kappa
density = peakHeightDensity( length(Dim)-1, kappa );
dz = 0.001;
z = min(locmaxZ):dz:max(locmaxZ); 
qcont = density(z);
%% %% Density of local maxima(discrete case)
rho_disc = exp(-1/2*1./(2*sigma_est.^2));

nz=length(z);
ddrho = 0.001;
drho=0:ddrho:0.999;
nrho=length(drho);
p=zeros(nz,nrho);
nnb=ones(nz,D-1)*2;
q=ddlm(z,ones(nz,D-1).*rho_disc,nnb,D-1);
qnormalize=q./(sum(q)*dz);

%% %% Fully connected DLM with sigma estimated from all
% locmaxZ3 = [];
% 
% sigma1 = discrete_covariance_ellipsoid(sigma_est, 3);
% niter = 1e5;
% n_lim = 1e6;
% 
% locmaxZ1 = [];
% nj = 10000;
% 
% for j = 1:nj
%     locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(3, sigma1, niter);
%     locmaxZ3 = [locmaxZ_sigma1 locmaxZ3];
%     if length(locmaxZ1) > n_lim
%         break;
%     end
% end

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

%% plot the 2D image with locmax highlighted
Z_loc = Z;
for i = 1:length(peakvals)
    Z_loc(peaklocs(1,i), peaklocs(2,i)) = (Z_loc(peaklocs(1,i), peaklocs(2,i))+5)*5;
end
imagesc(Z_loc);

%% generate multiplier fields and apply the method
imgs_2d = squeeze(imgs_cat(:,:,45,:));
multiplier = imgs_2d - mean(imgs_2d, 3);
mask_2d = logical(mask_2d);
base_fields = Field(multiplier, mask_2d);
new_field = multiplier_field( base_fields, nimgs);
new_field = Mask(new_field);
new_field = new_field.field;

subplot(3,1,1)
imagesc(base_fields.field(:,:,1))
subplot(3,1,2)
imagesc(mask_2d)
subplot(3,1,3)
imagesc(new_field(:,:,1))

%% Applying the white and gray matter
%mbs_loc = '/Users/tuo/Dropbox/Class HW and Notes/Rotation/2019_SP/code/real study/Gray.white_matter_masks';

fid1  = fopen('avg152T1_white.img'); %gray mask
data = fread(fid1);
GM_intensity = reshape(data, [91,109,91])/255;
%imagesc(GM_intensity(:,:,45));
colorbar

mask = GM_intensity > 0.5;
imagesc(mask(:,:,45))

% use the original images with the white matter mask
mask_2d = squeeze(mask(:,:,45));
imgs_white = zeros([Dim(1:2) nimgs]);
for I = 1:nimgs
    imgs_white(:,:,I) = squeeze(imgs_2d(:,:,I)).*mask_2d;
end

imgs_white = zero2nan(imgs_white);

% subtract the mean and standardize
mean_imgs_white =  mean(imgs_white, 3);
std_dev = std(imgs_white, 0, 3);
for I = 1:nimgs
    imgs_white(:,:,I) = (imgs_white(:,:,I) - mean_imgs_white)./std_dev;
end

% Get the localmax
locmaxZ2 = [];
locmaxZ = [];
B = 200;
nb = 40;
imgs_new = ones([Dim(1) Dim(2) B]);
load('inverse_rho_discrete.mat')

tic
for k = 1:B
    samp_ind = randsample(1:nimgs, nb);
    imgs_mask_2d_samp = imgs_white(:,:,samp_ind);
    Z = sqrt(nb)*mean(imgs_mask_2d_samp,3); %take the average and multiply by root-n
    Z_mask = Z(mask_2d == 1);
    Z_mask(isnan(Z_mask)) = 0;
    Z(mask_2d == 1) = Z_mask;
    [peaklocs, peakinds, peakvals] = lmindices(Z,"all",mask_2d);
    locmaxZ2 = [locmaxZ2;peakvals];
    
    % estimate the rho
    [fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est] = est_smooth(imgs_mask_2d_samp, mask_2d);
    
    sigma1 = discrete_covariance_ellipsoid(sigma_est, D-1);
    niter = 1e5;
    n_lim = 1e5;

    locmaxZ1 = [];
    nj = 10000;

    for j = 1:nj
        locmaxZ_sigma1 = simulateLocMax_Discrete_Stationary(D-1, sigma1, niter);
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
title("Beijing 4mm white matter mask")
ax = gca;

axis square

%% Histogram
ng = 1000;
histogram(locmaxZ2, 'Normalization','pdf', 'EdgeColor', '[0.9290 0.6940 0.1250]')
hold on
minz = min(locmaxZ);
maxz = max(locmaxZ);
pts = (minz-1):(maxz-minz+2)/ng:(maxz+1);
[g,xi] = ksdensity(locmaxZ, pts);
plot(xi, g, 'y')
plot(z,qcont,'g')
plot(z,qnormalize,'r')

legend('Peak of Image','Fully connected DLM','Continuous','Partially connected DLM',...
    'Location','northeast')

%% check the voxelwise distribution (narrower than normal?) check 198 and 40
mean_over_subjects_2d = mean_over_subjects(:,:,45);
std_dev_2d = std(imgs_2d, 0, 3);
i1 = 60;
i2 = 60;
histogram(imgs_2d(i1,i2,:), 'Normalization','pdf')
hold on
x = min(imgs_2d(i1,i2,:)):.01:max(imgs_2d(i1,i2,:));
y = normpdf(x,mean_over_subjects_2d(i1,i2),std_dev_2d(i1,i2));
plot(x, y)

ind = randsample(1:nimgs, 40);
imgs_2d_40 = imgs_2d(:,:,ind);
mean_40 = mean(imgs_2d_40, 3);
std_dev_40 = std(imgs_2d, 0, 3);
i1 = 40;
i2 = 40;
histogram(imgs_2d_40(i1,i2,:), 'Normalization','pdf')
hold on
x = [min(imgs_2d_40(i1,i2,:)):.01:max(imgs_2d_40(i1,i2,:))];
y = normpdf(x,mean_40(i1,i2),std_dev_40(i1,i2));
plot(x, y)

% try to use white noise with mask
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

imgs_white = zero2nan(imgs_white);

% subtract the mean and standardize
mean_imgs_white =  mean(imgs_white, 3);
std_dev = std(imgs_white, 0, 3);
for I = 1:nimgs
    imgs_white(:,:,I) = (imgs_white(:,:,I) - mean_imgs_white)./std_dev;
end

% Get the localmax
locmaxZ2 = [];
locmaxZ = [];
%B = 200;
B = nimgs;
nb = 100;
imgs_new = ones([Dim(1) Dim(2) B]);

tic
for k = 1:B
%     samp_ind = randsample(1:nimgs, nb);
%     imgs_mask_2d_samp = imgs_white(:,:,samp_ind);
%     Z = sqrt(nb)*mean(imgs_mask_2d_samp,3); %take the average and multiply by root-n
%       Z = imgs_white(:,:,k);
%       [peaklocs, peakinds, peakvals] = lmindices(Z,mask_2d);
    Z = imgs_white(:,:,k);
    [peaklocs, peakinds, peakvals] = lmindices(Z,'top',mask_2d);
%     Imax = imregionalmax(Z);
%     peakvals = Z(Imax);
    locmaxZ2 = [locmaxZ2;peakvals];
    
    % estimate the rho % 
    %[fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est] = est_smooth(imgs_mask_2d_samp, mask_2d);
    % estimate the rho
    Z_stand = Z/std(Z(~isnan(Z)));
    Z_mean = mean(Z(~isnan(Z)));
    rho_hor = (Z_stand(1:(Dim(1)-1),1:Dim(2))-Z_mean).*(Z_stand(2:Dim(1),1:Dim(2))-Z_mean);
    rho_vert = (Z_stand(1:Dim(1),1:(Dim(2)-1))-Z_mean).*(Z_stand(1:Dim(1),2:Dim(2))-Z_mean);
    %rho = [mean(rho_hor(~isnan(rho_hor))) mean(rho_vert(~isnan(rho_vert)))];
    %sigma_est = sqrt(1/4./-log(max(rho,10e-10))); % newton-raphson for
    %invrho fzero
    sigma_est = [0 0];
    
    %sigma1 = discrete_covariance_ellipsoid(sigma_est, 2);
    rho = discrete_covariance(0,D);
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

