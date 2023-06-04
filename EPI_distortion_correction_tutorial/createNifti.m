%**************************************************************************
%
%   Writes out data in nifti format.  If the data is complex then it is
%   written out as two files with the suffix _Magnitude & _Phase.
%
%   MFC 14.12.2012
%
%   MFC 01.05.2013  Oddly I had dat(:,:,:) = data so changed now to
%                   N.dat(:,:,:) = data
%
%   MFC 21.05.2013  Extra input parameter of spmMatrix so as to have
%                   orientational information.
%
%**************************************************************************

function createNifti(data, fn, spmMatrix)

% If a filename hasn't been passed in, where should I put it?
if nargin < 2
    [fn pn] = uiputfile('*.nii', 'Save As...');
    [~, fn] = fileparts(fn);
else
    [pn fn] = fileparts(fn);
    if isempty(pn)
        pn = pwd;
    end
        
end

% Is the data real or complex?
if isreal(data)
    
    dat         = file_array;
    dat.fname   = [pn filesep fn '.nii'];
    dat.dim     = size(data);
    dat.offset  = ceil(348/8)*8;
    
    if isfloat(data)
        dat.dtype   = 'FLOAT32';
        N.descrip = 'FLOAT data';
    else
        dat.dtype = 'UINT8';
        N.descrip = 'Non-float data';
    end
    
    % Create the NIFTI structure
    N = nifti;
    N.dat = dat;
    if exist('spmMatrix', 'var')
        N.mat = spmMatrix;
    end
    
    
    create(N); % Writes hdr
    N.dat(:,:,:,:) = data; % Writes data
    
else
    
    % General file_array fields:
    dat         = file_array;
    dat.dim     = size(data);
    dat.dtype   = 'FLOAT32';
    dat.offset  = ceil(348/8)*8;
    
    % Create the NIFTI structure for the magnitude:
    dat.fname   = [pn filesep fn '_Magnitude.nii'];
    N = nifti;
    N.dat = dat;
    N.descrip = 'Magnitude (from complex)';
    if exist('spmMatrix', 'var')
        N.mat = spmMatrix;
    end
    create(N); % Writes hdr info
    N.dat(:,:,:) = abs(data); % writes data to disk
    
    
    % Create the NIFTI structure for the phase:
    N = nifti;
    dat.fname = [pn filesep fn '_Phase.nii'];
    N.dat = dat;
    N.descrip = 'Phase (from complex)';
    create(N);
    N.dat(:,:,:) = angle(data);
    if exist('spmMatrix', 'var')
        N.mat = spmMatrix;
    end
    
end