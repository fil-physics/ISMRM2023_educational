function [data_corr , allAngle] = nyquist_ghostcorr(data, navs) 

%function [kdata_corr , allAngle] = nyquist_ghostcorr(kdata, kphasecor)
%
% Nadine N Graedel, n.graedel@ucl.ac.uk - WCHN/UCL 
% code provided for ISMRM 2023 educational lecture on EPI corrections
% last edit: 04/06/2023
%
% The following parameters need to be provided:
% - data: data to be corrected
% - navs: navigator data
%
% The data has the format: RO coils PE "segments", note that odd and even
% lines are stored in separate segments, the rest of the lines are zeros.
% Odd and even lines can be combined by adding over segments
% The navigators have the format: RO coils PE "averages" "segments",
% there are two "averages" because we acquire the one of the directions
% twice (three line navigator)
%
% "Segment 1" contains the data in the odd lines (zeros in even lines) =
% forward lines 
% "Segment 2" contains the data in the even lines (zeros in the odd lines)
% = reversed lines 
%  Two averages (as there are three phase correction lines) 

average_coils = true; % average correction over coils if desired
debug_plots   = true; % plot for example coil (coil 1)
point_by_point = false; % new point-by-point correction option (no linear assumption)
if point_by_point == true
    average_coils = false; 
    disp('WARNING: No average coil option for point by point phase correction');
end 

disp('Starting phase correction...');  
data_corr = zeros(size(data)); % output 
nCoil   = size(data,2);
nRO     = size(data,1);
allAngle = zeros(2,nCoil);
phs_diff_orig = zeros(nRO,nCoil);
    
% (1) get an estimated phase ramp for each coil
        for iCoil = 1:nCoil
            tmp_phasecor = navs(:,iCoil,:,:); 
            ft_phscor = ifftshift(fft(fftshift(tmp_phasecor,1),[],1),1);
            
            % (1.1) find number of points with signal greater or equal to at least half the maximum
            %signal 
            ft_phscor_ex = ft_phscor(:,1,1,1); % choose an example phae correction line 
            Nx = round(sum(abs(ft_phscor_ex)./ max(abs(ft_phscor_ex)) > 0.6)/2);
            %define range of points based on that 
            itv = nRO/2+1-Nx:nRO/2+1+Nx;
            xidx = -nRO/2:nRO/2-1;
        
            % (1.2) find which line is reflected (taking the mean is relevant for the
            %one where there are two lines...) 
            ft_phscor_fwd = squeeze(ft_phscor(:,1,1,1)); % take the odd line 
            ft_phscor_back = squeeze(mean(ft_phscor(:,1,:,2),3)); % take the mean of the even lines 

            % (1.3) find phase difference (complex division, then taking angle) 
            phs_diff_orig(:,iCoil)=angle( ft_phscor_fwd.* conj(ft_phscor_back) );
            magn_diff = abs(ft_phscor_fwd./ft_phscor_back);

            % (1.4) fit phase difference to straight line (polynomial of order 1) 
            angle_p = polyfit(xidx(itv).',phs_diff_orig(itv,iCoil),1);
            allAngle(:,iCoil) = angle_p;
            
            if debug_plots == true && iCoil == 1
                figure;
                tiledlayout(2,3); 
                nexttile; plot(xidx,abs(ft_phscor_fwd)); title('forward navigator'); legend('coil 1'); xlabel('samples'); ylabel('magnitude [au]');
                nexttile; plot(xidx,abs(ft_phscor_back)); title('backward navigators (averaged)'); xlabel('samples'); ylabel('magnitude [au]');
                nexttile; plot(xidx, magn_diff); title('magnitude ratio'); xlabel('samples'); ylabel('magnitude [au]');
                nexttile; plot(xidx,angle(ft_phscor_fwd)); title('phase'); xlabel('samples'); ylabel('phase [rad]');
                nexttile; plot(xidx,angle(ft_phscor_back)); title('phase'); xlabel('samples'); ylabel('phase [rad]');
                nexttile; plot(xidx,phs_diff_orig(:,iCoil)); title('phase difference'); hold all;  
                          plot(xidx,polyval(angle_p,xidx)); legend('phase diff','linear fit','Location','best'); xlabel('samples'); ylabel('phase [rad]');
                set(gcf,'color','w');
            end
        end
    
        % (2) Calculate the phase difference
        if (average_coils == true && point_by_point == false)
            disp('Averaging correction over coils...'); 
            % (2.1) In this case take the average over all coils (excluding
            % outlier coils)
            mean_angle = mean(allAngle,2); % average over coil dimension
            stdev_angle = std(allAngle,0,2); % std dev over coil dimension
     
            % Determine if outlier or not: 
            % If ||angle_coil - mean_angle|| > 2*std_angle over all coils, then this coil is
            % determined to be an outlier: 
            iNotOutlier = ~any(abs(allAngle - repmat(mean_angle,[1 nCoil]))>repmat(2*stdev_angle,[1 nCoil]),1);
     
            % for debugging: 
            tmp = ['# of outlier coils not used = ', num2str(nCoil - nnz(iNotOutlier))];
            disp(tmp); 

            %Take the mean accross all coils, which are not outliers: 
            mean_angle = mean(allAngle(:,iNotOutlier),2);
            phs_diff = polyval(mean_angle,xidx).'; % actually calculate the phase by evaluating the fit at the relevant locations (ones determined to have sufficient SNR)  
        
        else 
            disp('Correction determined and applied on a coil-by-coil basis...'); 
            phs_diff = zeros(length(xidx),nCoil);
            for iCoil = 1:nCoil
                phs_diff(:,iCoil) = polyval(allAngle(:,iCoil),xidx).';
            end
        end 
        
        % New point-by-point correction option (no linear assumption)
        if point_by_point == true
            phs_diff = phs_diff_orig;
        end 

        % (3) Apply the phase correction to the data
        for iCoil = 1:nCoil
             % Move each of forward and reverse lines to meet at the centre rather than
             % one line to meet the other (this is required in TURBINE sampling, otherwise
             % k-space will be asymmetric). For Cartesian 2D/3D EPI this only impacts the image phase

             img_fixtmp=ifftshift(fft(fftshift(data(:,iCoil,:,2),1),[],1),1); %reverse
             img_fixtmp2=ifftshift(fft(fftshift(data(:,iCoil,:,1),1),[],1),1); %forward

             sizeTemp = size(img_fixtmp);
             sizeTemp2 = size(img_fixtmp2);
             if average_coils == true 
                 phs_diffnew = repmat(phs_diff,[1,sizeTemp(2:end)]);
                 phs_diffnew2 = repmat(phs_diff,[1,sizeTemp2(2:end)]);
             else 
                 phs_diffnew = repmat(phs_diff(:,iCoil),[1,sizeTemp(2:end)]);
                 phs_diffnew2 = repmat(phs_diff(:,iCoil),[1,sizeTemp2(2:end)]);
             end 
             % apply phase shift in "image" space
             img_fixtmp = img_fixtmp.*exp(1i*phs_diffnew/2);
             img_fixtmp2 = img_fixtmp2.*exp(-1i*phs_diffnew2/2);
             
             data_corr(:,iCoil,:,2) = ifftshift(ifft(fftshift(img_fixtmp,1),[],1),1);  %reverse 
             data_corr(:,iCoil,:,1) = ifftshift(ifft(fftshift(img_fixtmp2,1),[],1),1); %forward
        end
              
disp('...completed phase correction.');  

end 


