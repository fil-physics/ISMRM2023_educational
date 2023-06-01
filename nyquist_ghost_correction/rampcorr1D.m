function kdata = rampcorr1D(kdata,t_ramp_up,t_flat_top,delay,ADC_dur)

% function kdata = rampcorr1D(kdata,t_ramp_up,t_flat_top,delay,ADC_dur)
%
% Nadine N Graedel, n.graedel@ucl.ac.uk - WCHN/UCL 
% code provided for ISMRM 2023 educational lecture on EPI corrections
% last edit: 31/05/2023
%
% The following parameters need to be provided:
% - kdata: readout must be first dimension
% - r_ramp_up: ramp up (and ramp down) time of the readout gradients
% - t_flat_top: flat top time of the readout gradients
% - delay: time at which the ADC is opened (relative to the start of the gradient)
% - ADC_dur: duration of sampling/ADC on

nRO           = size(kdata,1); 
t_ramp_down   = t_ramp_up;     % this script assumes symmetric trapezoid gradients
t             = linspace(delay, delay + ADC_dur, nRO); % time vector 
k_1D          = zeros(1, nRO);

%% (1) calculate k-space locations for ramp parameters given 
for iRO=1:nRO

% Explanation on the sampling calculation: 
% The goal is to calculate k_RO i.e. the area under the readout gradient. 
% We can assume G_max = 1 becasue this term cancels out with the normalization! 
% For simplicity subtract the intial triangle, during which there is no
% sampling (delay) in the end. 
        
    if t(iRO)<t_ramp_up % (a) sampling on the ramp up
        k_1D(iRO) = (0.5/t_ramp_up)*(t(iRO).^2);

    elseif t(iRO)>(t_ramp_up + t_flat_top)  % (b) sampling on the ramp down
        k_1D(iRO) = (0.5/t_ramp_up)*(t_ramp_up.^2) + (t(iRO)-t_ramp_up)-(0.5/t_ramp_down)*((t(iRO)-t_ramp_up-t_flat_top).^2);

    else % (c) sampling on the flat top
        k_1D(iRO) = (0.5/t_ramp_up)*(t_ramp_up.^2) + (t(iRO)-t_ramp_up); % full ramp up + amount of time spent on flat top so far... 
    end
end

%% (2) Regridding using a Sinc inperpolator 
OS_factor = 1; % nRO is determined from the data, therefore OS_factor = 1;
k1D_out = linspace(k_1D(1), k_1D(end), nRO/OS_factor)'; 
deltaK = k1D_out(2) - k1D_out(1);
    
density = diff(k_1D); density = [density density(1)]; 
grid_mat = sinc((repmat(k1D_out, 1, nRO)-repmat(k_1D,nRO/OS_factor, 1))/deltaK); 
grid_mat = repmat(density, nRO/OS_factor,1).*grid_mat;
grid_mat = grid_mat./(1e-12+repmat(sum(grid_mat,2),1, nRO));   
    
dims = size(kdata);
kdata = reshape(grid_mat*reshape(kdata, dims(1),[]),dims); % apply to all readout lines in this file 
        