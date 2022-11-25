%--------------------------------------------------------------------------
%% load 3d-gre: R9, 1mm
%--------------------------------------------------------------------------
% 60 Gb RAM is needed
clear
cd /project/3055010.04/RunningProjects/Harmonization/SharedBerkin/
addpath func/utils/
addpath func/utils_jose/

% sepia is used to access nii reading and writing functions
addpath('/home/common/matlab/sepia_1.0.0/')
sepia_addpath

% for HPC without access to matlab parallel computing where parallel jobs
% are submitted via fieldtrip
AtDondersHPC =0;
if AtDondersHPC ==1 ;
    addpath('/home/common/matlab/fieldtrip/qsub')
end
%% reconstruction of dataset without oversamplet

% first get all the sequencenames
MCNav_FullySampeld_data = 1;

Input.data_path = [pwd, '/AcquiredData/2022_10_28/'];
Input.DataFilename = ['1mm_clean/meas_MID01128_FID26878_pulseq_140_R9_1mm_mcNav.dat'];
Input.SeqFilename= ['SeqParamsUpdated_gre_3d_1mm_R9_mcNav.mat'];
Input.DataOutputFilename= ['gre_3d_R9_1mm_Rawdata'];


Options.Apply_B0_corr = 1 ;
Options.RemoveOs = 1 ;
Options.B0_corr_order = 0 ;

ReadOrderRemoveOsB0corrSeparateNav (Input,Options);



%% Creating Receive maps

Input2.data_path=Input.data_path;
Input2.DataInputFilename =Input.DataOutputFilename
Input2.FullKspaceDataVariable = 'NavFully'
Input2.DataOutputFilename ='Receive_R9'
CreateReceiveCoilsFromDataToData(Input2);

%% Recon data of the (no)B0 corrected Data
Nslices = 256
r2s = @ (ss,te) (real(1./ARLO(abs(ss),te(:))));
freq = @ (ss,dTE) (angle(mean(abs(ss(:,:,:,2:end)).*(ss(:,:,:,2:end)./ss(:,:,:,(1:end-1))),4))./(2*pi*(dTE)));
mkdir([Input.data_path,'recon'])

load([ Input.data_path ,Input.SeqFilename])

Input4.data_path = Input.data_path;
Input4.DataFilename= Input.DataOutputFilename;
Input4.DataName ='Data'
Input4.NavFilename= Input.DataOutputFilename;
Input4.NavName ='NavFully'
Input4.NavLabeln ='Nav_n'
Input4.ReceiveFilename= Input2.DataOutputFilename;
Input4.ReceiveName ='receive';

Options.lsqr_iter = 200;
Options.lsqr_tol = 1e-3;
Options.lambda = 1e-4;
Options.useNavdata = 1;


if AtDondersHPC == 1
    clear  Input_cell, clear Options_cell,
    SlicesPerJob = 8;
    for k = 1:Nslices/SlicesPerJob
        Input_cell{k} = Input4;
        Options.slices = ((k-1)*SlicesPerJob)+[1:SlicesPerJob];
        Options_cell{k} = Options;
    end
    
    Res_sum = qsubcellfun(@ReconstructTick, Input_cell,Options_cell,'memreq', 45 * 1024^3, 'timreq', 2000,'rerunable','yes','StopOnError',boolean(0));
    
    res = zeros([Nslices, size(Res_sum{1},2) ,size(Res_sum{1},3) ,size(Res_sum{1},4) ]);
    for k =1 :length(Res_sum)
        try
            res(Options_cell{k}.slices,:,:,:) = Res_sum{k};
        end
    end;
    
else
    delete(gcp('nocreate'))
    c = parcluster('local');        % build the 'local' cluster object
    
    total_cores = c.NumWorkers;
    parpool(ceil(total_cores * .5)) % use 50% of cores
    
    Options.slices =  1:NSlices             ;
    res = ReconstructTick(Input4, Options)
    delete(gcp('nocreate'))
    
end
save([ Input4.data_path , 'recon/',Input4.DataFilename,'_TikRecon'],'res' );
save_nii(make_nii(abs(res)),[ Input4.data_path , 'recon/',Input4.DataFilename,'_TikRecon.nii'] );


r2star = r2s(res,SeqParamsUpdated.TEff);

save([ Input4.data_path ,'recon/',Input4.DataFilename,'_TikRecon'],'r2star','-append' );


%% Recon data of the B0 corrected Data

clear Options
Options.lsqr_iter = 50;
Options.lsqr_tol = 1e-3;
Options.lambda = 1e-0;
Options.regul = 1;
Options.useNavdata = 0;
AtDondersHPC = 1;
if AtDondersHPC == 1
    clear  Input_cell, clear Options_cell,
    SlicesPerJob = 4;
    for k = 1:Nslices/SlicesPerJob
        Input_cell{k} = Input;
        Options.slices = ((k-1)*SlicesPerJob)+[1:SlicesPerJob];
        Options_cell{k} = Options;
    end

    Res_sum = qsubcellfun(@ReconstructTickJPM, Input_cell,Options_cell,'memreq', 30 * 1024^3, 'timreq', 2000,'rerunable','yes','StopOnError',boolean(0));

    res = zeros([Nslices, SeqParamsUpdated.Dims(2) ,SeqParamsUpdated.Dims(3) ,length(SeqParamsUpdated.TEff) ]);
    for k =1 :length(Res_sum)
        try
            res(Options_cell{k}.slices,:,:,:) =Res_sum{k};
        end
    end;

else
    delete(gcp('nocreate'))
    c = parcluster('local');        % build the 'local' cluster object

    total_cores = c.NumWorkers;
    parpool(ceil(total_cores * .5)) % use 50% of cores

    Options.slices =  1:Nslices             ;
    res = ReconstructTickJPM(Input, Options)
    delete(gcp('nocreate'))

end
% save([ Input.data_path ,Input.DataFilename,'_TikJPMRecon_NoPrior_vs2'],'res' );
save([Input4.data_path ,'recon/',Input4.DataFilename,'_TikJPMRecon_NoPrior_vs3'],'res' );

