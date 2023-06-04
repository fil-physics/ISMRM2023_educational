% Demo script Nyquist ghost correction
% Nadine N Graedel, n.graedel@ucl.ac.uk - WCHN/UCL 
% code provided for ISMRM 2023 educational lecture on EPI corrections
% last edit: 01/06/2023
%
%% load data 
% three data sets available for testing
% all data is unaccelerated 2D EPI
load('phantom_2DEPI_3T_1slc.mat');
%load('phantom_2DEPI_7T_1slc.mat');
%load('brain_2DEPI_7T_1slc.mat');

% The data has the format: RO coils PE "segments", note that odd and even
% lines are stored in separate segments, the rest of the lines are zeros.
% Odd and even lines can be combined by adding over segments
% The navigators have the format: RO coils PE "averages" "segments",
% there are two "averages" because we acquire the one of the directions
% twice (three line navigator)

%% 1D regridding to correct for ramp sampling
% comment out the two lines below if you want to see results with no ramp correction
data = rampcorr1D(data,ramp_up,flat_top,delay,ADCtime);
navs = rampcorr1D(navs,ramp_up,flat_top,delay,ADCtime);

%% Nyquist ghost correction
[data_corr , corr_factors_lin] = nyquist_ghostcorr(data, navs);
data_corr = sum(data_corr, 4); % combine odd and even lines 
data_nocorr = sum(data, 4);

%% FFT recon
im_corr = zeros(size(data_corr)); 
im_nocorr = zeros(size(data_nocorr)); 

for iCoil = 1:size(data_corr,2)
    im_corr(:,iCoil,:) = ifftshift(ifft(ifft(fftshift(data_corr(:,iCoil,:)),[],1),[],3));
    im_nocorr(:,iCoil,:) = ifftshift(ifft(ifft(fftshift(data_nocorr(:,iCoil,:)),[],1),[],3));
end

%% Sum of squares coil combination
im_corr_sos = squeeze(sqrt(sum(abs(im_corr).^2,2)));
im_nocorr_sos = squeeze(sqrt(sum(abs(im_nocorr).^2,2)));

%% Visualisation
% Figure 1 - plot k-space and image with/without correction for an example
iCoil = 1; % channels 1-32 

nRO = size(data_corr,1);
nPE = size(data_corr,3);
dPE = 15; % for zoom
dRO = 15; % for zoom

figure;
tiledlayout(2,3);
nexttile; imagesc(squeeze(abs(data_nocorr(:,iCoil,:)))); title('k-space single channel'); xlabel('PE'); ylabel('RO');
nexttile; imagesc(squeeze(abs(data_nocorr((nRO/2-dRO):(nRO/2+dRO),iCoil,(nPE/2-dPE):(nPE/2+dPE),:)))); title('k-space zoom'); xlabel('PE'); ylabel('RO');
nexttile; imagesc(squeeze(abs(im_nocorr(:,iCoil,:)))); title('image single channel'); xlabel('PE'); ylabel('RO');
nexttile; imagesc(squeeze(abs(data_corr(:,iCoil,:)))); title('k-space single channel corr'); xlabel('PE'); ylabel('RO');
nexttile; imagesc(squeeze(abs(data_corr((nRO/2-dRO):(nRO/2+dRO),iCoil,(nPE/2-dPE):(nPE/2+dPE),:)))); title('k-space zoom'); xlabel('PE'); ylabel('RO');
nexttile; imagesc(squeeze(abs(im_corr(:,iCoil,:)))); title('image single channel corr'); xlabel('PE'); ylabel('RO');
set(gcf,'position',[100,100,700,1000]) 
set(gcf,'color','w')

figure;
tiledlayout(1,2);
nexttile; imagesc(im_nocorr_sos(nRO/4+1:nRO*3/4,:)); title('sos image'); xlabel('PE'); ylabel('RO'); axis square;
nexttile; imagesc(im_corr_sos(nRO/4+1:nRO*3/4,:)); title('sos image with nyquist ghost correction'); xlabel('PE'); ylabel('RO'); axis square;
set(gcf,'position',[100,100,1500,600]); % position and size of the figure on the screen
set(gcf,'color','w');