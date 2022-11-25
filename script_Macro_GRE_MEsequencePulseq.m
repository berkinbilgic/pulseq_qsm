%--------------------------------------------------------------------------
%%  add path to relevant PulseqFiles
%--------------------------------------------------------------------------

CurrentDir = ['C:\Users\josmar\surfdrive\SharedBerkin']
cd(CurrentDir)
addpath([CurrentDir,'\pulseq-master\matlab'])
addpath([CurrentDir,'\func\utils_pulseq'])% some relvant functions to define k-space patterns

SaveDir =['compiled_seq_files/For2022Abstract/'];
mkdir(SaveDir)
%%
sys=mr.opts('maxGrad',24,'gradUnit','mT/m','riseTime',400e-6,...
    'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6,'B0',3);


%--------------------------------------------------------------------------
%% Protocol 1 = 3 mm iso, R=1, I went to 4mm to avoid having odd number and none integer matrix sizes
%--------------------------------------------------------------------------
DeltaTE = 4e-3;
TEmin = 3e-3;
TEmax = 32e-3;
res = 3;
SeqParams.TR  = 40e-3;                              % Repetition time in secs
SeqParams.TE  = [TEmin:DeltaTE:TEmax];              % Echo time in secs
SeqParams.Tread = 3.e-3;                            % Readouttime in secs
SeqParams.FlipAngle = 15;                           % flip angle in degrees
SeqParams.FOV = [256e-3] * [1 6/8 7/8] ;            % Field of View of output
SeqParams.Dims  = round( 256/res *[1 6/8 7/8]) ;    % Dimensions of output
SeqParams.useMCnavigator = 0 ;                      % in case a fully sampled centre of k-space
SeqParams.RyRzCaipi = [1 1 0] ;                     % [Ry Rz Caipi] should define the downsampling pattern
SeqParams.Ndummy = round(5 / SeqParams.TR)  ;       % 5 secs of Dummys
SeqParams.rfphasecycle = 50;                        % we are using the 50 degree used in the siemens implementatoin
SeqParams.readoutOversampling  =2;                  % readout Oversampling siemens uses Factor

SeqParams.OutputFile = [SaveDir, 'gre_3d_',num2str(res),'mm_R',num2str(prod(SeqParams.RyRzCaipi(1:2))),'.seq']     % Output File

tic
[seq, SeqParamsUpdated] = writeGradientEcho3D_ME_FreeSpacing_CAIPI (SeqParams, [])
toc

Display_MainAtributes(seq,SeqParamsUpdated)

%--------------------------------------------------------------------------
%% Protcol 2 : 3mm acquisition where each echo was downsampled by a factor of 9
%--------------------------------------------------------------------------


DeltaTE = 4e-3;
TEmin = 3e-3;
TEmax = 32e-3;
res = 3;
SeqParams.TR  = 40e-3;                              % Repetition time in secs
SeqParams.TE  = [TEmin:DeltaTE:TEmax];              % Echo time in secs
SeqParams.Tread = 2.5e-3;                           % Readouttime in secs
SeqParams.FlipAngle = 15;                           % flip angle in degrees
SeqParams.FOV = [256e-3] * [1 6/8 7/8] ;            % Field of View of output
SeqParams.Dims  = round( 256/res *[1 6/8 7/8]);     % Dimensions of output
SeqParams.useMCnavigator = 0                        % in case a fully sampled centre of k-space
SeqParams.RyRzCaipi = [1 9 3] ;                     % [Ry Rz Caipi] should define the downsampling pattern
SeqParams.Ndummy = round(5 / SeqParams.TR)  ;       % 5 secs of Dummys
SeqParams.rfphasecycle = 50;                        % we are using the 50 degree used in the siemens implementatoin
SeqParams.readoutOversampling  =2;                  % readout Oversampling siemens uses Factor
SeqParams.OutputFile = [SaveDir, 'gre_3d_',num2str(res),'mm_R',num2str(prod(SeqParams.RyRzCaipi(1:2))),'.seq']     % Output File

tic
[seq, SeqParamsUpdated] = writeGradientEcho3D_ME_FreeSpacing_CAIPI (SeqParams, [])
toc

Display_MainAtributes(seq,SeqParamsUpdated)
%%

%--------------------------------------------------------------------------
%% Protcol 2 : 3mm acquisition where each echo was downsampled by a factor of 9
%--------------------------------------------------------------------------

DeltaTE = 4e-3;
TEmin = 3e-3;
TEmax = 32e-3;
res = 3;
SeqParams.TR  = 40e-3;                             % Repetition time in secs
SeqParams.TE  = [TEmin:DeltaTE:TEmax];                                    % Echo time in secs
SeqParams.Tread = 2.5e-3;                            % Readouttime in secs
SeqParams.FlipAngle = 15;                              % flip angle in degrees
SeqParams.FOV = [256e-3] * [1 6/8 7/8] ;             % Field of View of output
SeqParams.Dims  = round( 256/res *[1 6/8 7/8])  ;         % Dimensions of output
SeqParams.useMCnavigator = 1;  % in case a fully sampled centre of k-space
SeqParams.RyRzCaipi = [1 9 3] ;     % [Ry Rz Caipi] should define the downsampling pattern
SeqParams.Ndummy = round(5 / SeqParams.TR)  ;   % 5 secs of Dummys
SeqParams.rfphasecycle = 50    ;                 % we are using the 50 degree used in the siemens implementatoin
SeqParams.readoutOversampling  =2;                  % readout Oversampling siemens uses Factor
SeqParams.OutputFile = [SaveDir, 'gre_3d_',num2str(res),'mm_R',num2str(prod(SeqParams.RyRzCaipi(1:2))),'_mcNav.seq']     % Output File

tic
[seq, SeqParamsUpdated] = writeGradientEcho3D_ME_FreeSpacing_CAIPI (SeqParams, sys)
toc

Display_MainAtributes(seq,SeqParamsUpdated)
%%

%--------------------------------------------------------------------------
%% Protcol 3 : 3mm acquisition where each echo was downsampled by a factor of 9 + echo shifting used to
%--------------------------------------------------------------------------

DeltaTE = 4e-3;
TEmin = 3e-3;
TEmax = 32e-3;
res = 3;
GoldenAngleFraction = 0.38; % I will use this as a way to unalias successive pairs of echos

SeqParams.TR  = 40e-3;                               % Repetition time in secs
SeqParams.Tread = 2.5e-3;                            % Readouttime in secs
SeqParams.TE  = [TEmin:DeltaTE:TEmax];                                    % Echo time in secs
for k=  3 : 2 :length(SeqParams.TE)
    SeqParams.TE(k:end) = SeqParams.TE(k:end) + round(DeltaTE * GoldenAngleFraction*1e5)/1e5;
end
SeqParams.FlipAngle = 15;                           % flip angle in degrees
SeqParams.FOV = [256e-3] * [1 6/8 7/8] ;            % Field of View of output
SeqParams.Dims  = round( 256/res *[1 6/8 7/8]) ;    % Dimensions of output
SeqParams.useMCnavigator = 0;                       % in case a fully sampled centre of k-space
SeqParams.RyRzCaipi = [1 9 3] ;                     % [Ry Rz Caipi] should define the downsampling pattern
SeqParams.Ndummy = round(5 / SeqParams.TR)  ;       % 5 secs of Dummys
SeqParams.rfphasecycle = 50 ;                       % we are using the 50 degree used in the siemens implementatoin
SeqParams.readoutOversampling  =2;                  % readout Oversampling siemens uses Factor
SeqParams.OutputFile = [SaveDir, 'gre_3d_',num2str(res),'mm_R',num2str(prod(SeqParams.RyRzCaipi(1:2))),'_varES.seq']     % Output File

tic
[seq, SeqParamsUpdated] = writeGradientEcho3D_ME_FreeSpacing_CAIPI (SeqParams, sys)
toc

Display_MainAtributes(seq,SeqParamsUpdated)

%--------------------------------------------------------------------------
%% Protcol 4 : 1mm acquisition where each echo was downsampled by a factor of 9
%--------------------------------------------------------------------------

DeltaTE = 5e-3;
TEmin = 3e-3;
TEmax = 32e-3;
res = 1;
SeqParams.TR  = 40e-3;                                  % Repetition time in secs
SeqParams.Tread = 3e-3;                               % Readouttime in secs
SeqParams.TE  = [TEmin:DeltaTE:TEmax];                                    % Echo time in secs
SeqParams.FlipAngle = 15;                               % flip angle in degrees
SeqParams.FOV = [256e-3] * [1 6/8 7/8] ;                % Field of View of output
SeqParams.Dims  = round( 256/res *[1 6/8 7/8]) ;        % Dimensions of output
SeqParams.useMCnavigator = 1;                           % in case a fully sampled centre of k-space
SeqParams.RyRzCaipi = [1 9 3] ;                         % [Ry Rz Caipi] should define the downsampling pattern
SeqParams.Ndummy = round(5 / SeqParams.TR)  ;           % 5 secs of Dummys
SeqParams.rfphasecycle = 50 ;                           % we are using the 50 degree used in the siemens implementatoin
SeqParams.readoutOversampling  = 2;                      % readout Oversampling siemens uses Factor
SeqParams.OutputFile = [SaveDir, 'gre_3d_',num2str(res),'mm_R',num2str(prod(SeqParams.RyRzCaipi(1:2))),'.seq']     % Output File

tic
[seq, SeqParamsUpdated] = writeGradientEcho3D_ME_FreeSpacing_CAIPI (SeqParams, sys)
toc

Display_MainAtributes(seq,SeqParamsUpdated)

%%

SeqParams.OutputFileMat = [SaveDir, 'gre_3d_',num2str(res),'mm_R',num2str(prod(SeqParams.RyRzCaipi(1:2))),'.mat']

load (SeqParams.OutputFileMat);

% mc_order = reshape(SeqParamsEff.navcount,[max(SeqParamsEff.navcount) , length(SeqParamsEff.navcount)/(max(SeqParamsEff.navcount))]);
mc_order = reshape(SeqParamsEff.navcount,[ length(SeqParamsEff.navcount)/(max(SeqParamsEff.navcount)), max(SeqParamsEff.navcount) ]);

mc = mc_order(round(end/2),:);
Ny = SeqParamsEff.DimsNav(2);
Nz = SeqParamsEff.DimsNav(3);
nTE = length(SeqParamsEff.TEff)

Nav_n =zeros([Ny,Nz,nTE]);
Acq = 0;
iShotData = 0;
iShotNav = 0;
iReadOutNav = 0;
for iShot = 1:length(SeqParamsEff.labelData0_Nav1)
    
    if SeqParamsEff.labelData0_Nav1(iShot) == 0
        iShotData = iShotData + 1;
        
    else
        iShotNav = iShotNav + 1;
        [PE1, PE2, TE] = ind2sub([Ny,Nz,nTE],find(KspaceOrderNavKyKzt==iShotNav));
        
        for t = 1:nTE
            
            Acq = Acq + 1;
            iReadOutNav = iReadOutNav + 1;
            Nav_n(PE1(t), PE2(t),t) = SeqParamsEff.navcount (iReadOutNav);
        end
        
    end
end
Nav_mask =zeros([Ny,Nz,10]);
k= 0;
for nav = mc(1:10)
    k =k+1
    for  t = 1:nTE
        Nav_mask(:,:,k)=Nav_mask(:,:,k)+ t*double((Nav_n(:,:,t)==nav));
    end;
end


for k = 1:9
    subplot(3,3,k)
    imab(crop(Nav_mask(:,:,k),[9 9]*2))
    title(['k - space mask MC nav ',num2str(k) ])
    
end;
colormap(gray)

