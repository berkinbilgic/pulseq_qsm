%--------------------------------------------------------------------------
%% load 3d-gre: R1, 4mm
%--------------------------------------------------------------------------

data_path = '/autofs/space/marduk_001/users/berkin/2022_08_17_bay1_pulseq_gre_invivo/';

    
% filename1 = 'meas_MID00310_FID16272_pulseq_140_R1_3mm_cor_R__L.dat';        % mtx 85x64x75
% filename1 = 'meas_MID00312_FID16274_pulseq_140_R9_3mm_cor_R__L.dat';        % mtx 85x72x81

% filename = 'meas_MID00313_FID16275_pulseq_140_R9_3mm_cor_R__L_mcNav.dat';
% filename = 'meas_MID00314_FID16276_pulseq_140_R9_3mm_cor_R__L_varES.dat';

% filename = 'meas_MID00316_FID16278_pulseq_140_R4_1mm_cor_R__L.dat';
filename = 'meas_MID00318_FID16280_pulseq_140_R9_1mm_cor_R__L.dat';

tic
    dt = mapVBVD([data_path, filename]);
    res = dt{end}.image.unsorted();
toc



%--------------------------------------------------------------------------
%% load 3d-gre
%--------------------------------------------------------------------------

save_path = '/autofs/cluster/berkin/berkin/Matlab_Code_New/pulseq/HarmonizationMEGREPulseq_v0/compiled_seq_files/20220421/';

% filename = 'gre_3d_3mm_R1';
% filename = 'gre_3d_3mm_R9';
% filename = 'gre_3d_3mm_R9_varES';
% filename = 'gre_3d_3mm_R9_mcNav';

filename = 'gre_3d_1mm_R9';
% filename = 'gre_3d_1mm_R4';


load([save_path, 'SeqParamsUpdated_', filename, '.mat'])


num_chan = size(res,2);
num_echoes = length(SeqParamsUpdated.TE);
Dims = SeqParamsUpdated.Dims;

Data = zeros([Dims,num_chan,num_echoes]);
NavFully = zeros([Dims(1:3),num_chan,num_echoes]);


Ny = Dims(2);
Nz = Dims(3);
nTE = num_echoes;

KspaceMaskKyKzt = SeqParamsUpdated.KspaceMaskKyKzt;

if isfield(SeqParamsUpdated, 'KspaceOrderNavKyKzt')
    disp('mcNav acquired')
    KspaceOrderNavKyKzt = SeqParamsUpdated.KspaceOrderNavKyKzt;
else
    disp('no mcNav')
end

iShotData = 0;
iShotNav = 0;
Acq = 0;

for iShot = 1:length(SeqParamsUpdated.labelData0_Nav1)
    if ~(mod(iShot,100))
        disp(['parsing shot ', num2str(iShot), ' / ', num2str(length(SeqParamsUpdated.labelData0_Nav1))])
    end
    
    if SeqParamsUpdated.labelData0_Nav1(iShot) == 0
        iShotData = iShotData + 1;
        [PE1, PE2, TE] = ind2sub([Ny,Nz,nTE],find(KspaceMaskKyKzt==iShotData));

        for t = 1:num_echoes
            Acq = Acq + 1;
            Data(:, PE1(t), PE2(t),:,t) = res(:,:,Acq);
        end
        
    else
        iShotNav = iShotNav + 1;
        [PE1, PE2, TE] = ind2sub([Ny,Nz,nTE],find(KspaceOrderNavKyKzt==iShotNav));
        
        for t = 1:num_echoes
            Acq = Acq + 1;

            NavFully(:, PE1(t), PE2(t),:,t) = res(:,:,Acq);
        end

    end
end


%--------------------------------------------------------------------------
%% 3 mm iso case:
% R9 FOV -> 511.2, 216.5, 243.6     (3mm iso, 85x72x81 mtx)
% R1 FOV -> 511.2, 192.4, 225.5     (3mm iso, 85x64x65 mtx)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% 1mm iso case:
% R4 FOV -> 512, 192, 224           (1mm iso, 512x192x224 mtx)
% R9 FOV -> 512, 198, 225           (1mm iso, 512x198x225 mtx)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% crop and zero-pad mcnav from R=9 case
%--------------------------------------------------------------------------

% reference data: eco1 from R9 mcnav acquisition
load([data_path, 'img_ref_eco1_R9_mcnav'])

% img_ref = img_nav(:,:,:,:,1);

s_img_ref = s(img_ref(:,:,:,1));

% pad to 1mm resolution
img_ref_pad = ifft3call(padarray(fft3call(img_ref), s_img_ref));

% match fov
img_ref_pad = padarray(img_ref_pad, [1,0,0], 'post');

% R4 1mm case
% img_ref_pad = img_ref_pad(:,1+12:end-12,1+10:end-9,:);

% R9 1mm case
img_ref_pad = img_ref_pad(:,1+9:end-9,1+9:end-9,:);


%--------------------------------------------------------------------------
%% remove 2x oversampling
%--------------------------------------------------------------------------

img = ifft3call(Data);
img = img(1+end/4:3*end/4,:,:,:,:);


if isfield(SeqParamsUpdated, 'KspaceOrderNavKyKzt')
    img_nav = ifft3call(NavFully);
    img_nav = img_nav(1+end/4:3*end/4,:,:,:,:);
end


for t = 1:num_echoes
    imagesc3d2(rsos(img(:,:,:,:,t),4), s(img)/2, t, [-90,180,180], [0.,1e-4]), setGcf(.5)
end

R_factor = zeros(num_echoes,1);

for t = 1:num_echoes
    tmp = Data(:,:,:,1,t);
    R_factor(t) = 1/mean(tmp(:)~=0);
end

disp(R_factor)


%--------------------------------------------------------------------------
%% esprit  
%--------------------------------------------------------------------------

addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/TOOLBOXES/Espirit_matlab_only_toolbox'))

num_acs = 24;
kernel_size = [6,6];
eigen_thresh = 0.7;             % threshold for coil sens mask size

delete(gcp('nocreate'))
c = parcluster('local');        % build the 'local' cluster object

total_cores = c.NumWorkers;     % 48 cores for marduk
parpool(ceil(total_cores * .75))


% img_R1_pad = img_nav(:,:,:,:,1);
img_R1_pad = img_ref_pad(:,:,:,:,1);

receive = zeross(s(img_R1_pad(:,:,:,:,1)));

tic    
parfor slc_select = 1:s(img_R1_pad,1)
    disp(num2str(slc_select))
    
    kspace_slice = fft2call(sq(img_R1_pad(slc_select,:,:,:,1)));
        
    [maps, weights] = ecalib_soft( kspace_slice, num_acs, kernel_size, eigen_thresh );

    receive(slc_select,:,:,:) = dot_mult(maps, weights >= eigen_thresh); 
end
toc

delete(gcp('nocreate'))

% save([data_path, 'receive_3mm_R9_mcnav.mat'], 'receive', '-v7.3')
% save([data_path, 'receive_1mm_R4_mcnav.mat'], 'receive', '-v7.3')
save([data_path, 'receive_1mm_R9_mcnav.mat'], 'receive', '-v7.3')



%--------------------------------------------------------------------------
%% load coil maps
%--------------------------------------------------------------------------

% load([data_path, 'receive_3mm_R9_mcnav.mat'])
load([data_path, 'receive_1mm_R4_mcnav.mat'])


%--------------------------------------------------------------------------
%% sense recon
%--------------------------------------------------------------------------

addpath '/autofs/cluster/berkin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox'


lsqr_iter = 200;
lsqr_tol = 1e-3;

img_res = zeross(s(sq(img(:,:,:,1,:))));
m2d = sq(Data(1,:,:,:,:)) ~=0;

mosaic(sq(m2d(:,:,1,:)),2,4,56,'',[0,1]),colormap parula
mosaic(mean(m2d(:,:,1,:),4),1,1,57,'',[0,1]),colormap parula

    
delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; % 48 cores for marduk
parpool(ceil(total_cores/2))

tic
parfor slc_select = 1:s(img,1)
    disp(num2str(slc_select))
    
    sens = sq(receive(slc_select,:,:,:));
    kspace_slc = fft2call(sq(img(slc_select,:,:,:,:)));


    param = [];
    param.N = Dims(2:3);
    param.num_chan = num_chan;
    param.lambda = 1e-4;        % L2 reg

    param.sens = sens;

    Res = zeross([param.N, num_echoes]);

    for t = 1:num_echoes
        kspace_coils = kspace_slc(:,:,:,t);
        param.m2d = m2d(:,:,:,t);

        res = lsqr(@apply_sense_tikc, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  
        
        Res(:,:,t) = reshape(res, param.N);
    end

    img_res(slc_select,:,:,:) = Res;

%     mosaic(sq(img_res(slc_select,:,:,:)), 2, ceil(num_echoes/2), 16, num2str(slc_select), [0,1e-3], 90)
end
toc

delete(gcp('nocreate'))

% save([data_path, 'img_res_3mm_R9.mat'], 'img_res', '-v7.3')
% save([data_path, 'img_res_3mm_R9_varES.mat'], 'img_res', '-v7.3')
% save([data_path, 'img_res_3mm_R9_mcNav.mat'], 'img_res', '-v7.3')

% save([data_path, 'img_res_1mm_R4.mat'], 'img_res', '-v7.3')
% save([data_path, 'img_res_1mm_R9.mat'], 'img_res', '-v7.3')


for t = 1:num_echoes
    imagesc3d2(img_res(:,:,:,t), s(img_res)/2, t, [-90,180,-180], [0.,6e-4]), setGcf(.5)
end

mosaic(padarray(sq(m2d(1:32,1:32,1,:)),[4,4]),2,3,11,'',[0,1]),setGcf(.5), colormap parula
mosaic(mean(sq(m2d(1:32,1:32,1,:)),3),1,1,12,'',[0,1]),setGcf(.5), colormap parula


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

img_nav_combo = zeross(size(img_res));

for t = 1:num_echoes
    img_nav_combo(:,:,:,t) = coil_combine(img_nav(:,:,:,:,t), receive);
    
    imagesc3d2(img_nav_combo(:,:,:,t), s(img_res)/2, t, [90,90,-90], [0.,4e-3]), setGcf(.5)
end



%--------------------------------------------------------------------------
%% load recons
%--------------------------------------------------------------------------

load([data_path, 'img_res_3mm_R9_varES.mat']);  img_R9_res_varES = img_res;
load([data_path, 'img_res_3mm_R9.mat']);        img_R9_res = img_res;


for t = 1:2:num_echoes
    imagesc3d2(img_R9_res_varES(:,:,:,t), s(img_res)/2, t, [90,90,-90], [0.,4e-3]), setGcf(.5)
    imagesc3d2(img_R9_res(:,:,:,t), s(img_res)/2, t+10, [90,90,-90], [0.,4e-3]), setGcf(.5)
end



