function CreateReceiveCoilsFromDataToData (Input)
% 
% Input.data_path,  
% Input.DataInputFilename
% Input.DataOutputFilename
% Input.FullKspaceDataVariable = 'NavFully' 


%--------------------------------------------------------------------------
%% esprit: parfor
%--------------------------------------------------------------------------

% addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/TOOLBOXES/Espirit_matlab_only_toolbox'))
 

NavFully =load([Input.data_path,  Input.DataInputFilename] , Input.FullKspaceDataVariable)
NavFully = single(NavFully.(Input.FullKspaceDataVariable));

num_acs = 24;
kernel_size = [6,6];
eigen_thresh = 0.7;             % threshold for coil sens mask size

% delete(gcp('nocreate'))
% c = parcluster('local');        % build the 'local' cluster object
% 
% total_cores = c.NumWorkers;
% parpool(ceil(total_cores * .5)) % use 50% of cores

img_nav = ifft3call(NavFully);
clear NavFully 

img_R1_pad = img_nav(:,:,:,:,1);
clear img_nav
receive = zeross(s(img_R1_pad(:,:,:,:,1)));

tic

  for slc_select = 1:s(img_R1_pad,1)
%  for slc_select = 20

    disp(num2str(slc_select))

    kspace_slice = fft2call(sq(img_R1_pad(slc_select,:,:,:,1)));

    [maps, weights] = ecalib_soft( kspace_slice, num_acs, kernel_size, eigen_thresh );

    receive(slc_select,:,:,:) = dot_mult(maps, weights >= eigen_thresh);

end

toc

% delete(gcp('nocreate'))
try
display(' trying appending')
save([Input.data_path,  Input.DataOutputFilename],  'receive','-append','-v7.3' )
display(' appending done')
catch
display(' trying a new file')
save([Input.data_path,  Input.DataOutputFilename],  'receive','-v7.3' )
display(' saving done')
end
