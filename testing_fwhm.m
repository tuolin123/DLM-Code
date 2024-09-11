contrast = 'WM';
global dataloc
subjectids = filesindir([dataloc,'HCP\HCPContrasts\', contrast, '\']);
data = zeros([91,109,91,80]);
for I = 1:length(subjectids)
    %modul(I,10)
    subject_file_loc = [dataloc, 'HCP\HCPContrasts\WM\',...
                                    subjectids{I}, '\WM\Level2\cope11.nii.gz'];
    data(:,:,:,I) = imgload(subject_file_loc);
end

MNImask = imgload('MNImask');

%%
est_smooth(data,0,MNImask)

%%
[rho, sigma] = estsd(data - mean(data), mask, 3, 0.5);
sigma2FWHM(sigma) 

%%
sb_dir = statbrainz_maindir;
fid1  = fopen([sb_dir, 'BrainImages\Volume\avg152T1_gray.img']);
GMdata = fread(fid1);
GM_intensity = reshape(GMdata, [91,109,91])/255;
GM_mask = GM_intensity > 0.5;
imagesc(GM_mask(:,:,50))

%%
sb_dir = statbrainz_maindir;
fid1  = fopen([sb_dir, 'BrainImages\Volume\avg152T1_white.img']);
WMdata = fread(fid1);
WM_intensity = reshape(WMdata, [91,109,91])/255;
WM_mask = WM_intensity > 0.5;
imagesc(WM_mask(:,:,50))

%%
est_smooth(data,0,WM_mask)

%%
[rho, sigma] = estsd(data - mean(data), GM_mask, 3, 0.5);
sigma2FWHM(sigma)

%%
Dim = [91,109,91];
nsubj = 80;
FWHM = 3;
noise = noisegen(Dim, nsubj, FWHM);
D = 3;
initial = 0.5;
% mask = zeros(Dim); mask(30:60,40:80,30:60) = 1; mask = logical(mask);
mask = GM_mask;
%noise = noise.*GM_mask;
%noise = noise - mean(noise);
[rho_est, sigma_est] = estsd(noise, GM_mask, D, initial);
sigma2FWHM(sigma_est) 

%%
est_smooth(noise, 0, GM_mask)