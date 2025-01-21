%% read data 
hcp = spm_read_vols(spm_vol('cope11.nii.gz'));
dim = size(hcp);
D = length(dim) - 1;

% not masking
mask = spm_read_vols(spm_vol('MNImask.nii'));
mask = dilate_mask(mask, -1);

FWHM = 0;
isotropic = 0;
gaussian = 0; % set this to 1 if gaussianize
niter = 1e5;

% estimate the sd for smoothing kernel
[fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est] = est_smooth(hcp, mask);
fwhm_kernel = sigma2FWHM(sigma_est);

% identify peaks and calculate peak p values by MCDLM
tic
[peaklocs, peakinds, peakvals, peakpvals] = MCDLM(hcp, mask, D, FWHM, isotropic, gaussian, niter);
toc

ind_min = find(peakpvals == min(peakpvals));

pval = peakpvals(~isnan(peakpvals));
pmin = pval(find(pval == min(pval)));

sum(pval<0.05)+ind_min-1

% continuous method
kappa   = 1; 
density = peakHeightDensity( D, kappa );
dz = 0.0001;
z = -6:dz:14;
qcont = density(z);
pval_cont = 1-cumsum(qcont(2:end,:).*dz);
peak_t2normal = -norminv(tcdf(-peakvals,80)); % Gaussianization
n_inf = sum(peak_t2normal == inf);
pval_hcp_cont = zeros(length(peak_t2normal),1);
for i = (n_inf+1):length(peak_t2normal)
    dif = z-peak_t2normal(i);
    ind = find(abs(dif) == min(abs(dif)));
    pval_hcp_cont(i) = pval_cont(ind);
end
sum(pval_hcp_cont<0.05)

%plot figures
hcp_tfield = convfield_t(hcp, 0);

brainDim1 = squeeze(hcp_tfield.field(32,:,:));
brainDim1 = brainDim1.';
brainDim2 = squeeze(hcp_tfield.field(:,68,:));
brainDim2 = brainDim2.';
brainDim3 = squeeze(hcp_tfield.field(:,:,65));
brainDim3 = brainDim3.';
t = tiledlayout(1,3);
nexttile
imagesc(brainDim1);
set(gca,'YDir','normal')
xlim([11 101]);
axis off
nexttile
imagesc(brainDim2);
set(gca,'YDir','normal')
axis off
nexttile
imagesc(brainDim3);
set(gca,'YDir','normal')
ylim([11 101]);

brainDim1 = squeeze(hcp_tfield.field(20,:,:));
brainDim1 = brainDim1.';
brainDim2 = squeeze(hcp_tfield.field(:,43,:));
brainDim2 = brainDim2.';
brainDim3 = squeeze(hcp_tfield.field(:,:,62));
brainDim3 = brainDim3.';
t = tiledlayout(1,3);
nexttile
imagesc(brainDim1);
set(gca,'YDir','normal')
xlim([11 101]);
nexttile
imagesc(brainDim2);
set(gca,'YDir','normal')
nexttile
imagesc(brainDim3);
set(gca,'YDir','normal')
ylim([11 101]);

%Adjusting form multiple comparisons - BH
ind_min = find(peakpvals == min(peakpvals));
peakpvals(1:(ind_min-1)) = min(peakpvals);
peakpvals(isnan(peakpvals)) = max(peakpvals);
alpha = 0.05;
[ rejection_ind, nrejections, sig_locs ] = fdrBH(peakpvals, alpha);

%MCDLM with Gaussianization - need to load the p-value first
ind_min = find(peakpvals == min(peakpvals));
peakpvals(1:(ind_min-1)) = min(peakpvals);
peakpvals(isnan(peakpvals)) = max(peakpvals);
alpha = 0.05;
[ rejection_ind1, nrejections1, sig_locs1 ] = fdrBH(peakpvals, alpha);

%continuous
[ rejection_ind2, nrejections2, sig_locs2 ] = fdrBH(pval_hcp_cont, alpha);


