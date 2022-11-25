%% 20180216 KC: Fixed fv doesn't work for fsl version>5.0.9 in the cluster
function fv(data,varargin)
data = double(data);
% function fv(data,varargin)
% varargin
% 'voxdims',[] default is [1 1 1],
% 'sepvols',0/1

 
[voxdims, bSepVols, isfsleyes] = process_options(varargin,'voxdims',[],'sepvols',0,'fsleyes',false);

if strcmp(computer,'MACI64')
    setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
end
 
if strcmp(computer,'GLNXA64')
    tempbase = '~/temp';
    % try to read the fslversion from FSL root directory
    [~,fsldir] = unix('echo $FSLDIR');
    ind = strfind(fsldir,'/');
    fslversion = str2double(fsldir(ind(end)+5:end));
    fslMajorRelease = str2double(fsldir(ind(end)+1));
    if (fslMajorRelease==5 && fslversion > 9) || fslMajorRelease > 5
        if isfsleyes
            fslview = 'fsleyes';
        else
            fslview = 'fslview_deprecated';
        end
    else
        fslview = 'fslview';
    end
% else
%     tempbase = '~/temp';
end
if ~isdir(tempbase),
   mkdir (tempbase);
end;
        unix(['rm ~/temp/temp* ']);

tempFile = [tempbase '/temp_fslview'];
tempFile_phs = [tempbase '/temp_fslview_phs'];

 
if (bSepVols)
    % if using separate volumes - take the absolute of the input data...
    if ~isreal(data(1))
        data = abs(data);
    end
    nVols = size(data,4);
    sysTxt = 'fslview';
    for iV = 1:nVols
        thisTxt = [tempFile '_' num2str(iV,'%.3d') '.nii'];        
        save_nii(make_nii(data(:,:,:,iV),voxdims),thisTxt);
        sysTxt = [sysTxt ' ' thisTxt];
    end
    system([sysTxt ' &']);    

    
else
    if isreal(data)
        save_nii(make_nii(data,voxdims),[tempFile '.nii']);
%         unix(['fslview ' tempFile ' &']);
        unix([fslview ' ' tempFile ' &']);
    else
        save_nii(make_nii(abs(data),voxdims),tempFile);
        save_nii(make_nii(angle(data),voxdims),[tempFile_phs '.nii']);
%         unix(['fslview ' tempFile_phs ' ' tempFile ' &']);
        unix([fslview ' ' tempFile_phs ' ' tempFile ' &']);
        
    end
end

