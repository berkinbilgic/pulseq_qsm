
function   [ seq, SeqParamsEff]  = writeGradientEcho3D_ME_CAIPI (SeqParams,sys)

%%
%[ seq, SeqParams]  = writeGradientEcho3D_ME_CAIPI (SeqParams,sys)
% Seqparams is the update seqparams in case TR and TE had to be updated, it
% also has information on the
%
% SeqParams.TR              Repetition time in secs
% SeqParams.TE              Echo time in secs - it does not have to be
% defined explicitely, you can choose to only define Tread
% SeqParams.Tread           Readouttime in secs
% SeqParams.FlipAngle       flip angle in degrees
% SeqParams.FOV             Field of View of output
% SeqParams.Dims            Dimensions of output
% SeqParams.OutputFile      Output File
% SeqParams.useMCnavigator  in case a fully sampled centre of k-space
% SeqParams.RyRzCaipi       [Ry Rz Caipi] should define the downsampling pattern
% SeqParams.FreqPe1Pe2      string 'xyz' or 'zyx' and so on X is LR, Y AP, Z HF
% SeqParams.Ndummy          default is ~5secs of dummies
% SeqParams.rfphasecycle    the rf phase cycle phase desired (in degrees) 
%%
% clear
SeqParamsEff = SeqParams;

if isfield(SeqParams,'TR');
    TR = SeqParams.TR;
else
    TR = 40e-3; % in secs
end

if isfield(SeqParams,'FlipAngle');
    FlipAngle = SeqParams.FlipAngle;
else
    FlipAngle = 8; % in secs
end

if isfield(SeqParams,'FOV');
    FOV = SeqParams.FOV;
else
    FOV=[256e-3] * [1 6/8 7/8];     % Define FOV in m; this is something tht would cover the whole head in the saggital direction
end


if isfield(SeqParams,'Dims');
    dims = SeqParams.Dims;
else
    dims=round( 64 *[1 6/8 7/8]);
end

if isfield(SeqParams,'TE');
    TE = SeqParams.TE;
    nTE = length(TE);
    defineEchoTimes = true;
    
else
    defineEchoTimes = false;
    TE = [];
end


if isfield(SeqParams,'Tread')
    Tread = SeqParams.Tread;
    defineReadoutTime = true;
else
    defineReadoutTime = false;
    Tread = [];
end

if isfield(SeqParams,'RyRzCaipi')
    Ry = SeqParams.RyRzCaipi(1);
    Rz = SeqParams.RyRzCaipi(2);
    CaipiShift = SeqParams.RyRzCaipi(3);
else
    Ry = 1;
    Rz = 2;
    CaipiShift = 1;
end

if isfield(SeqParams,'useMCnavigator')
    useMCnavigator = SeqParams.useMCnavigator;
else
    useMCnavigator = 0;
end

if isfield(SeqParams,'OutputFile')
    OutputFile = SeqParams.OutputFile
else
    OutputFile = 'gre3d.seq';
end

% SimParams.FreqPe1Pe2      string 'xyz' or 'zyx' and so on X is LR
if isfield(SeqParams,'FreqPe1Pe2')
    FreqPe1Pe2 = SeqParams.FreqPe1Pe2;
else
    FreqPe1Pe2 = 'xyz';
end
if isfield(SeqParams,'Ndummy')
    Ndummy = SeqParams.Ndummy ;
else
    Ndummy = ceil(5 / TR); % if it is not defined it will take at least 5 seconds
end

if isfield(SeqParams,'rfphasecycle')
    rfphasecycle = SeqParams.rfphasecycle ;
else
    rfphasecycle = 117/2; % if it is not defined it will used the 117/2 that was initially on the pulseq code
end



% Nx = dims(1); Ny = dims(2) ,Nz = dims(3);
res = FOV/dims;
%% alternative options to decide on the number of echo times and ADC times
% defineEchoTimes = false;
% defineReadoutTime = true;
% useMCnavigator = 1; % only partially introduced the use of a navigator for image correction

if and(defineEchoTimes,isempty(Tread))
    Tread = [TE(1)/2]; %initialization is half the first echo time
end

%%
%% Define down sampling pattern, the given pattern defines the new effective FOV Nx, Ny, Nz
if useMCnavigator
    res_nav = 8e-3; % resolution associated with the navigator  16mm for convenience
    dimsNav =ceil(dims./(res_nav/res)); % this defines the navigator region
end

% Ry = 1; Rz = 12; CaipiShift = 5; % acceleration of each echo
[pat] = makeKernelTight(Ry, Rz, CaipiShift);

dimspat = [1 size(pat)];            %
% figure(1)
% imab (pat)
dimsNew = ceil(dims./dimspat).*dimspat
if useMCnavigator
    dimsNavNew = ceil(dimsNav./dimspat).*dimspat;
    SeqParamsEff.DimsNav = dimsNavNew;
    
end

% update FOV if needed
FOV = dimsNew .*res;
Nx = dimsNew(1); Ny = dimsNew(2) ;Nz = dimsNew(3);

SeqParamsEff.Dims = dimsNew;
SeqParamsEff.FOV = FOV;


%%

% define system properties
if isempty(sys)
    riseTime=400e-6;                %
    sys=mr.opts('maxGrad',24,'gradUnit','mT/m','riseTime',riseTime,...
        'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6,'B0',3);
end;
% what the systems returns as maximum gradient is gamma * maxGrad
seq=mr.Sequence(sys);           % Create a new sequence object


% Create non-selective pulse
[rf, rfDelay] = mr.makeBlockPulse(FlipAngle*pi/180,sys,'Duration',0.2e-3);

% Define other gradients and ADC events
deltak=1./FOV;


gx = mr.makeTrapezoid(FreqPe1Pe2(1),'FlatArea',Nx*deltak(1),'FlatTime',Tread,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
adc.dwell = round(adc.dwell /seq.adcRasterTime) *seq.adcRasterTime; % ensures the timings work...
% gx = mr.makeTrapezoid(FreqPe1Pe2(1),'FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys);
% adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);

%  this assumes it's a blip
% Tpre = ceil(2*sqrt(deltak*Nx/2/sys.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time and assume a blip will be used
% gxPre = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'Area',-gx.area/2,'Duration',Tpre(1));

% computes the fastest possible prephased and rephaser, this will later be
% relaxed to minimize eddy currents
gxPre = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'Area',-gx.area/2,'amplitude',sys.maxGrad);
Tpre = mr.calcDuration(gxPre);
gxRep = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'Area',-gx.area,'amplitude',sys.maxGrad);
Trep = mr.calcDuration(gxRep);

% gxSpoil = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'Area',gx.area,'Duration',TRep);
% JPM changed to be used as a frequency navigator, it has to cross kspace
% center again, to it will have twie the are of gx
gxSpoil = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'Area',gx.area*2,'Duration',Tread/2); % arbitrarilly ask the spoiler to have half the duration of the readout
% for some reason this wasy to prescribe using the flat are does not work
% gxSpoil = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'FlatArea',-gx.flatArea*2,'Duration',round(Tread/2/seq.gradRasterTime)*seq.gradRasterTime); % note that my navigator might not be fully centred in k-space... but it will nevertheless have the same resolution as the normal readouts, and, because we will be only looking at differences over time that should be ok..

TSpoil = mr.calcDuration(gxSpoil);

% adc_1dnav = mr.makeAdc(Nx,'Duration',(gxSpoil.flatTime)/2,'Delay',gxSpoil.riseTime,'system',sys);
adc_nav = mr.makeAdc(Nx,'Duration',(gx.amplitude/gxSpoil.amplitude)*gx.flatTime,'Delay',gxSpoil.riseTime); 
adc_nav.dwell = round(adc_nav.dwell /seq.adcRasterTime) *seq.adcRasterTime; % ensures the timings work...


% check how many echo times fit in the TR
% sequence consists of:
% rf delay gxPre   nTE [ gx gxPre] GSpoil)

% TimeToFillWithEchotimes =floor((TR - mr.calcDuration(rf) - mr.calcDuration(gxPre)  ...
%     - mr.calcDuration(gxSpoil))/seq.gradRasterTime)*seq.gradRasterTime;
TimeToFillWithEchotimes =floor((TR - rf.delay - mr.calcDuration(rf) - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gxSpoil))/seq.gradRasterTime)*seq.gradRasterTime;
if defineReadoutTime
    % checks how many echo times would fit in the available space
    nTEmax = floor(TimeToFillWithEchotimes /(mr.calcDuration(gxRep) + mr.calcDuration(gx)))
    
    if defineEchoTimes
        if nTEmax<nTE
            display('the prescribed TEs were not possible given the readout time')
            display('the TEs will be adapted using the defined Tread times')
            nTE = nTEmax;
            TE = mr.calcDuration(rf)/2 + mr.calcDuration(gxPre) + mr.calcDuration(gx)/2 + [0:(nTE-1)] * (mr.calcDuration(gx) + Trep );
        else
            display('The TEs perscribed were possible to create')
        end;
        Trep = floor((min(diff(TE)) - mr.calcDuration(gx))/seq.gradRasterTime)*seq.gradRasterTime;
    else
        % TRrep (time to fill with echo times / number of echo times minus the readout time - it should be the time to apply the refocusing)
        nTE = nTEmax;
        Trep = floor((TimeToFillWithEchotimes/nTE - mr.calcDuration(gx))/seq.gradRasterTime)*seq.gradRasterTime;
        TE = mr.calcDuration(rf)/2 + mr.calcDuration(gxPre) + mr.calcDuration(gx)/2 + [0:(nTE-1)] * (mr.calcDuration(gx) + Trep );
    end;
end
%

% Calculate timing
delayTE = floor((TE(1) - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
% delayTR = ceil((TR - mr.calcDuration(rf) - mr.calcDuration(gxPre) ...
%     - (mr.calcDuration(gx) + Trep ) *nTE - delayTE - mr.calcDuration(gxSpoil))/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = ceil((TR - mr.calcDuration(rfDelay) - mr.calcDuration(gxPre) ...
    - (mr.calcDuration(gx) + Trep ) *nTE - delayTE - mr.calcDuration(gxSpoil))/seq.gradRasterTime)*seq.gradRasterTime;

TReff = ((mr.calcDuration(rfDelay) + mr.calcDuration(gxPre) ...
    + (mr.calcDuration(gx) + Trep ) *nTE + delayTE + delayTR + mr.calcDuration(gxSpoil)));

if delayTE*delayTR<0 % checks if ones of this is negative, a possibility when explicitely defining TEs
    if delayTE<0
        TE = TE - delayTE + seq.gradRasterTime; % the TEs had to be increased
        %         delayTR = delayTR + delayTE - seq.gradRasterTime;
        delayTR = round((delayTR + delayTE - seq.gradRasterTime)/seq.gradRasterTime)*seq.gradRasterTime;
        delayTE = seq.gradRasterTime;
        display('the TEs will be adapted using the defined Tread times')
        
    end;
    if delayTR<0
        TE = TE + delayTR - seq.gradRasterTime; % the TEs had to be decreased
        %         delayTE = delayTR + delayTE - seq.gradRasterTime;
        delayTE  = round((delayTR + delayTE - seq.gradRasterTime)/seq.gradRasterTime)*seq.gradRasterTime;
        delayTR = seq.gradRasterTime;
        display('the TEs will be adapted using the defined Tread times and TR')
        
    end;
end

dTE = mr.makeDelay(delayTE);
dTR = mr.makeDelay(delayTR);


SeqParamsEff.res = res;
SeqParamsEff.TR = TReff;
SeqParamsEff.TE = TE;
SeqParamsEff.Tread = Tread;
SeqParamsEff.TRdelay = delayTR;
SeqParamsEff.TEdelay = delayTE;
SeqParamsEff.Trf = mr.calcDuration(rf);
SeqParamsEff.Trep = Trep;
SeqParamsEff.Tpre = Tpre;
SeqParamsEff.Ndummy = Ndummy
%% compute sampling patterns and trajectories
% this part is for the full k-space data
KspaceMaskKyKzt = kspacepattern_filling_order ([Nx,Ny,Nz,nTE],pat,[floor(Ry/2) round(CaipiShift/2)] ); % the shift across echos is a fraction of the caipishift
labelData0_Nav1 = zeros ([1 max(KspaceMaskKyKzt(:))]); % 0 given to Data, 1 given to MCnavdata
SeqParamsEff.KspaceMaskKyKzt = KspaceMaskKyKzt;

% this appraoch works well in the case where pat was defined by putting all
% the acceleration along the z direction, and thus large caipi jumps might
% be used there

if useMCnavigator
    %         keyboard
    
    [KspaceOrderNavKyKzt , navcount] = kspacepattern_filling_navorder([dimsNavNew,nTE],pat);
    % HERE this still has to be introduced in a spread manner in the normal
    % acquisition, the k-space also has to be zero padded before further work
    %     creates a vector from 0 - 1
    imaging_MCNavs = linspace(0,1,max(KspaceOrderNavKyKzt(:)));
    imaging_Data = linspace(0,1,max(KspaceMaskKyKzt(:)));
    
    %      KspaceOrderNavKyKzt = KspaceOrderNavKyKzt + max(KspaceMaskKyKzt(:)); % this can currently be done in this way because the nav data is fully sampled
    KspaceOrderNavKyKzt = zpad(KspaceOrderNavKyKzt,size(KspaceMaskKyKzt));
    % finding correct order, how to spread the navigators through the imaging process
    [~, shot_ordering] = sort(cat(2,imaging_Data,imaging_MCNavs));
    labelData0_Nav1 = double(shot_ordering > max(KspaceMaskKyKzt(:))); % 0 given to Data, 1 given to MCnavdata
    SeqParamsEff.KspaceOrderNavKyKzt = KspaceOrderNavKyKzt;
    SeqParamsEff.navcount = navcount;
end;

SeqParamsEff.labelData0_Nav1 = labelData0_Nav1;
%%
%
% keyboard

% Make trapezoids for inner loop to save computation
clear gyPre gyReph;
% defines the area of various blips
areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);
% defines k-space coordinate defined from there index1 = kmin
jumpY = ((0:Ny-1)-Ny/2);
jumpZ = ((0:Nz-1)-Nz/2);

gxRep = mr.makeTrapezoid(FreqPe1Pe2(1),sys,'Area',-gx.area,'Duration',Trep);
for iY=1:Ny
    gyPre(iY) = mr.makeTrapezoid(FreqPe1Pe2(2),'Area',areaY(iY),'Duration',Tpre);
    gyReph(iY) = mr.makeTrapezoid(FreqPe1Pe2(2),'Area',-areaY(iY),'Duration',Trep); %
end
for iZ=1:Nz
    gzPre(iZ) = mr.makeTrapezoid(FreqPe1Pe2(3),'Area',areaZ(iZ),'Duration',Tpre);
    gzReph(iZ) = mr.makeTrapezoid(FreqPe1Pe2(3),'Area',-areaZ(iZ),'Duration',Trep);
end
for iY=1:Ny
    gyPre(iY).id = seq.registerGradEvent(gyPre(iY));
    gyReph(iY).id = seq.registerGradEvent(gyReph(iY));
end
for iZ=1:Nz
    gzPre(iZ).id = seq.registerGradEvent(gzPre(iZ));
    gzReph(iZ).id = seq.registerGradEvent(gzReph(iZ));
end


% preregister constant objects to accelerate computations
% this is not necessary, but accelerates the sequence creation by up to a factor of 2
% there is one more place in the second loop
gxPre.id=seq.registerGradEvent(gxPre);
gx.id=seq.registerGradEvent(gx);
gxSpoil.id=seq.registerGradEvent(gxSpoil);
gxRep.id = seq.registerGradEvent(gxRep);


%adc.id=seq.registerAdcEvent(adc);
%dTE.id=seq.registerDelayEvent(dTE);
%dTR.id=seq.registerDelayEvent(dTR);
%rfDelay.id=seq.registerDelayEvent(rfDelay);
[~, rf.shapeIDs]=seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes

% Drive magnetization to the steady state
for iY=1:Ndummy
    % RF
    rf.phaseOffset = mod(rfphasecycle*(iY^2+iY+2)*pi/180,2*pi);
    seq.addBlock(rf,rfDelay);
    % Gradients
    seq.addBlock(gxPre,gyPre(round(Ny/2)),gzPre(round(Nz/2)));
    seq.addBlock(dTE);
    %     seq.addBlock(gx);
    %     seq.addBlock(gyReph(floor(Ny/2)),gxSpoil);
    
    for iE = 1:nTE
        seq.addBlock(gx);
        if iE ~= nTE
            seq.addBlock(gxRep)
        else
            seq.addBlock(gyReph(round(Ny/2)),gzReph(round(Nz/2)))
            
        end
    end;
    seq.addBlock(gxSpoil);
    
    seq.addBlock(dTR);
end

iShotData = 0;
iShotNav = 0;

for iShot=1:length(labelData0_Nav1)
    
    
    if labelData0_Nav1(iShot) == 0
        iShotData = iShotData + 1;
        [ijumpY ,ijumpZ, pois] = ind2sub([Ny,Nz,nTE],find(KspaceMaskKyKzt==iShotData));
        %         length(ijumpY);
        
    else
        iShotNav = iShotNav + 1;
        [ijumpY ,ijumpZ, pois] = ind2sub([Ny,Nz,nTE],find(KspaceOrderNavKyKzt==iShotNav));
        
        %         length(ijumpY);
    end
    
    %     gzPre = mr.makeTrapezoid(FreqPe1Pe2(z),'Area',areaZ(iZ),'Duration',Tpre);
    %     gzReph = mr.makeTrapezoid(FreqPe1Pe2(z),'Area',-areaZ(iZ),'Duration',Tpre);
    %     % optional pre-registration for acceleration
    %     gzPre.id = seq.registerGradEvent(gzPre);
    %     gzReph.id = seq.registerGradEvent(gzReph);
    % RF spoiling
    rf.phaseOffset = mod(rfphasecycle*(iShot^2+iShot+2)*pi/180/2,2*pi); % correction
    adc.phaseOffset = rf.phaseOffset;
    adc_nav.phaseOffset  = rf.phaseOffset;
    % Excitation
    seq.addBlock(rf,rfDelay);
    
    % Encoding
    %     seq.addBlock(gxPre,gyPre(ijumpY(1)),gzPre(ijumpZ(1)));
    PrePhaserBlockContents={gxPre,gyPre(ijumpY(1)),gzPre(ijumpZ(1))}; % here we demonstrate the technique to combine variable counter-dependent content into the same block
    PrePhaserBlockContents= { PrePhaserBlockContents{:} , mr.makeLabel('SET','ECO', 0), mr.makeLabel('SET','NAV',false) }; % set the nav counters
    PrePhaserBlockContents= { PrePhaserBlockContents{:}, mr.makeLabel('SET','LIN', ijumpY(1)-1), mr.makeLabel('SET','PAR', ijumpZ(1)-1) }; % set the lin and partitions
    if labelData0_Nav1(iShot) == 0
        PrePhaserBlockContents= { PrePhaserBlockContents{:} , mr.makeLabel('SET','SET', 0)} ;
    else
        PrePhaserBlockContents= { PrePhaserBlockContents{:} , mr.makeLabel('SET','SET', navcount((iShotNav-1) * nTE + 1 ))};
        
    end;
    
    seq.addBlock(PrePhaserBlockContents{:});
    
    
    seq.addBlock(dTE);
    for iE = 1:nTE
        seq.addBlock(gx,adc);
        %         seq.addBlock(gx);
        if iE ~= nTE
            
            djumpY = -ijumpY(iE+1)+ijumpY(iE); % blip needed in kspace to go for the next step, polarity is unexpected because rephase is done with the negative area
            djumpZ = -ijumpZ(iE+1)+ijumpZ(iE);
            
            %           seq.addBlock( mr.makeLabel('INC','ECO', 1),...
            %           gxRep,gyReph(find(jumpY==djumpY)),gzReph(find(jumpZ==djumpZ)))
            %           ; % should move to the next k-space line and change the echo
            %           time
            RephaserBlockContents = {gxRep,gyReph(find(jumpY==djumpY)),gzReph(find(jumpZ==djumpZ))};
            RephaserBlockContents = {RephaserBlockContents{:}, ...
                                        mr.makeLabel('SET','ECO', iE), ...
                                        mr.makeLabel('SET','LIN', ijumpY(iE+1)-1), ...
                                        mr.makeLabel('SET','PAR', ijumpZ(iE+1)-1)};
            if labelData0_Nav1(iShot) ~= 0
                RephaserBlockContents= { RephaserBlockContents{:} , mr.makeLabel('SET','SET', navcount((iShotNav-1) * nTE +iE + 1 ))};
                
            end;
            seq.addBlock(RephaserBlockContents{:});
        else
            %         seq.addBlock( gyReph(ijumpY(iE)), gzReph(ijumpZ(iE)) ) ;
            seq.addBlock(mr.makeLabel('SET','NAV',true),gxRep, gyReph(ijumpY(iE)), gzReph(ijumpZ(iE)) ) ;
        end
    end;
    seq.addBlock(gxSpoil,adc_nav);
    %      seq.addBlock(gxSpoil);
    seq.addBlock(dTR);
end
SeqParamsEff.adc = adc;


fprintf('Sequence ready\n');

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;
%

SeqParamsEff.error_report = error_report;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Visualise sequence and output for execution
seq.plot('TimeRange',[Ndummy+1 Ndummy+3]*TR)

seq.setDefinition('FOV', FOV);
seq.setDefinition('Name', 'gre3d_me_nav');
seq.write(OutputFile,false);

% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
if Nx<=64
    tic;
    [kfa,ta,kf,t,texcite]=seq.calculateKspacePP(); % ta is the timing of measurements and kfa the k-space coordinate at those times
    toc
    figure;
    subplot(111)
    plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
    
    plot3(ta(:),kfa(2,:),kfa(3,:),'r.');
    
    figure;
    %     set(gcf,'Position',[2.6697    0.1043    0.5600    0.4200]*1000)
    
    subplot(121)
    
    plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
    
    subplot(122)
    plot3(ta(:),kfa(2,:)/deltak(2),kfa(3,:)/deltak(3),'r.');
    xlabel('Acquisition Time')
    ylabel('Ky number')
    zlabel('Kz number')
end


%% create a smoothly rotating plot
if Nx<=16
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
    kabsmax=max(abs(kf)')';
    kxyabsmax=max(kabsmax(1:2));
    kxyzabsmax=max(kabsmax);
    %axis([-kxyabsmax kxyabsmax -kxyabsmax kxyabsmax -kabsmax(3) kabsmax(3)])
    axis([-kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax])
    [caz,cel] = view;
    for caz_add=0:1:359
        view(caz+caz_add,cel);
        drawnow;
    end
end

%% create a smoothly rotating plot (rotated to read along z)
if Nx<=16
    figure;plot3(kf(2,:),-kf(3,:),kf(1,:));
    hold on;plot3(kfa(2,:),-kfa(3,:),kfa(1,:),'r.');
    set(gca,'visible','off'); % hide axes
    set(gca, 'CameraViewAngle',get(gca, 'CameraViewAngle')); % freeze the view
    kabsmax=max(abs(kf)')';
    kxyabsmax=max(kabsmax(1:2));
    kxyzabsmax=max(kabsmax);
    %axis([-kxyabsmax kxyabsmax -kxyabsmax kxyabsmax -kabsmax(3) kabsmax(3)])
    %axis([-kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax])
    s1=1.2;
    axis([ -kabsmax(2)*s1 kabsmax(2)*s1  -kabsmax(2)*s1 kabsmax(2)*s1 min(kf(1,:)) kabsmax(1)]);
    [caz,cel] = view;
    folder='kspace3d';
    mkdir(folder);
    for caz_add=0:1:359
        view(caz+caz_add,cel);
        drawnow;
        print( '-r100', '-dpng', [folder '/frame_' num2str(caz_add,'%03d') '.png']);
        % use convert frame_???.png -gravity center -crop 300x300+0+0 +repage -delay 0.1 -loop 0 kspace_gre3d.gif
        % to create a GIF movie
    end
end
