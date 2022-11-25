function ReadOrderRemoveOsB0corrSeparateNav (Input,Options);

% Input.data_path = [pwd, '/AcquiredData/2022_09_invivo_data/'];
% Input.DataFilename = ['meas_MID00313_FID16275_pulseq_140_R9_3mm_cor_R__L_mcNav.dat'];
% Input.SeqFilename= ['gre_3d_3mm_R9_mcNav'];
% Input.DataOutputFilename= ['gre_3d_3mm_R9_mcNav'];
% Options.Apply_B0_corr =1 or 0
% Options.RemoveOs;
% Options.B0_corr_order = 0 (default) 1 first order ;

tic
dt = mapVBVD([Input.data_path, Input.DataFilename]);
%     res = dt{end}.image.unsorted();
toc

if isfield(Options,'B0_corr_order')
    fitfirstorder = Options.B0_corr_order;
else
    fitfirstorder = 0;
end
%--------------------------------------------------------------------------
%% load 3d-gre
%--------------------------------------------------------------------------


load([Input.data_path, Input.SeqFilename])



num_echoes = length(SeqParamsUpdated.TE);
Dims = SeqParamsUpdated.Dims;

%  keyboard

%% extracting the respiration navigation


if Options.Apply_B0_corr == 1 ;

    res = dt{end}.phasecor.unsorted();


    % taking data to real space and removing oversampling
    res_ft = fftshift(fft(fftshift(res, 1), [], 1), 1);
    res_ft = res_ft ./ mean(res_ft, 3).* abs(mean(res_ft, 3));
    if Options.RemoveOs
        res_ft = res_ft(1+end/4:3*end/4,:,:,:,:);
    end;

    % creating regression model

%     fitfirstorder = 0;

    if fitfirstorder
        model=double(zeros(size(res_ft,1),2));
        model(:,1)=1;
        model(:,2)=(1:size(res_ft,1))-size(res_ft,1)/2;%x -siemens
        %fit first order term
        b  = zeros ( 2, size(res_ft,2), size(res_ft,3));
    else
        model=double(zeros(size(res_ft,1),1));
        model(:,1)=1;
        %fit first order term
        b  = zeros ( 1, size(res_ft,2), size(res_ft,3));
    end;
    for shot = 1:size(res_ft,3)
        for chan = 1:size(res_ft,2)
            R = angle( res_ft(:, chan, shot));
            w = abs( res_ft(:, chan, shot)).^1; % We are emphasizing regions with signal
            w( w < max(w)/10) = 0; % and ignoring where signal is under a 5th of maximum signal
            b(:,chan,shot) = pinv(bsxfun(@times,model,w))*(R.*w);
            %         b(:,chan,shot) = pinv(bsxfun(@times,model,w))*(R.*w);

        end;
    end;
    %

    AverageWindowIter = 20
    B0_smooth = b ;
    for k = 1 :AverageWindowIter
        B0_smooth(:,:,2:end-1) = 1/3 * (B0_smooth(:,:,1:end-2) + B0_smooth(:,:,2:end-1) + B0_smooth(:,:,3:end) );
    end;
    % checking the fit works well
    chan = 16;
    figure
    subplot(221)
    plot(squeeze(angle(res_ft(:,:,shot))) + [1:32])
    hold on
    plot(model*b(:,chan,shot)+ [1:32],'-.')
    title('fit of phase change of profile for the 32 channels')
    xlabel('position')

    subplot(222)
    hold off
    plot(squeeze(b(1,chan,:)))
    hold on
    plot(squeeze(B0_smooth(1,chan,:)),'r')
    
    plot(find(SeqParamsUpdated.labelData0_Nav1==0),squeeze(b(1,chan,SeqParamsUpdated.labelData0_Nav1==0)),'g')
    try
    plot(find(SeqParamsUpdated.labelData0_Nav1==1),squeeze(b(1,chan,SeqParamsUpdated.labelData0_Nav1==1)),'k')
    end
    xlabel('shot')
    title('timecourse of zero order fit')
    axis tight
    if fitfirstorder
        subplot(224)
        plot(squeeze(b(2,chan,:)))
        hold on
        plot(squeeze(B0_smooth(2,chan,:)),'r')
        plot(find(SeqParamsUpdated.labelData0_Nav1==0),squeeze(B0_smooth(2,chan,SeqParamsUpdated.labelData0_Nav1==0)),'g')
        xlabel('shot')
        title('timecourse of first order fit')
        axis tight
    end
    clear res_ft

    figure
    subplot(313)
    plot(squeeze(mean(b(1,chan,:),2)))
    hold on
    plot(squeeze(mean(B0_smooth(1,chan,:),2)),'r')
    xlabel('N excitation')
    ylabel('zero order correction')
    title('all data and small smooth')
    subplot(312)
    plot(find(SeqParamsUpdated.labelData0_Nav1==0),squeeze(mean(b(1,chan,SeqParamsUpdated.labelData0_Nav1==0),2)),'g')
    ylabel('zero order correction')

    title('Imaging data alone')

    subplot(311)
    try
    plot(find(SeqParamsUpdated.labelData0_Nav1==1),squeeze(mean(b(1,chan,SeqParamsUpdated.labelData0_Nav1==1),2)),'k')
    ylabel('zero order correction')
    title('Nav data alone')
    end
    if fitfirstorder
        figure

        subplot(313)
        plot(squeeze(mean(b(2,chan,:),2)))
        ylabel('first order correction')

        xlabel('N excitation')
        hold on
        plot(squeeze(mean(B0_smooth(2,chan,:),2)),'r')
        title('all data and small smooth')
        subplot(312)
        plot(find(SeqParamsUpdated.labelData0_Nav1==0),squeeze(mean(b(2,chan,SeqParamsUpdated.labelData0_Nav1==0),2)),'g')
        ylabel('first order correction')

        title('Imaging data alone')
try
        subplot(311)
        plot(find(SeqParamsUpdated.labelData0_Nav1==1),squeeze(mean(b(2,chan,SeqParamsUpdated.labelData0_Nav1==1),2)),'k')
        ylabel('first order correction')
        title('Nav data alone')
end
    end;
end;
%%
%%
% keyboard
% removes oversampling at the start
res = dt{end}.image.unsorted();
num_chan = size(res,2);
res = fftshift(fft(fftshift(res, 1), [], 1), 1);
if Options.RemoveOs
    res = res(1+end/4:3*end/4,:,:,:,:);
end;
if Options.Apply_B0_corr == 1
    % estimation of Navigator echo time
    TE_nav = 0.5*((SeqParamsUpdated.TR -SeqParamsUpdated.TRdelay)...
        + SeqParamsUpdated.TEff(end));
    for shot = 1:size(B0_smooth,3)
        for chan = 1:num_chan
            phasecorr = -model*B0_smooth(:,chan,shot);
            for te = 1:num_echoes
                res(:,chan,(shot-1)*num_echoes +te) = res(:,chan,(shot-1)*num_echoes +te) ...
                    .*exp(1i * phasecorr *SeqParamsUpdated.TEff(te)/TE_nav);

            end;
        end;
    end;
end;
% goes back to k-space
res = ifftshift(ifft(ifftshift(res, 1), [], 1), 1);

Dims(1) = Dims(1)/2;

Data = single(zeros([Dims,num_chan,num_echoes]));
NavFully = single(zeros([Dims(1:3),num_chan,num_echoes]));
Nav_n = zeros([1, Dims(2:3),num_chan,num_echoes]);
% this matrix allows to extract the MCnavs from the fully encoded NavFully
% by doing
% mcnav_n = 1
% MCnav = NavFully .* double(Nav_n == mcnav_n);



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
iReadOutNav = 0; % added to be able to correctly atribute readouts to MC navigators

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
            Data(:, PE1(t), PE2(t),:,t) = single(res(:,:,Acq));
        end

    else
        iShotNav = iShotNav + 1;
        [PE1, PE2, TE] = ind2sub([Ny,Nz,nTE],find(KspaceOrderNavKyKzt==iShotNav));

        for t = 1:num_echoes

            Acq = Acq + 1;
            iReadOutNav = iReadOutNav + 1;
            NavFully(:, PE1(t), PE2(t),:,t) = single(res(:,:,Acq));
            Nav_n(1, PE1(t), PE2(t),:,t) = SeqParamsUpdated.navcount(iReadOutNav);
        end

    end
end
m2d = sq(Data(1,:,:,:,:)) ~=0;


%% code to extract one single mcNav
% first mask the KspaceMaskKyKzt
%readouts associated with mcnav_n
% mcnav_n = 1
% MCnav = NavFully .* double(Nav_n == mcnav_n);

clear KspaceOrderNavKyKzt


%%
clear dt
clear res
clear b
% clear SeqParamsUpdated

%% --------------------------------------------------------------------------

% img = ifft3call(Data);
% % img = img(1+end/4:3*end/4,:,:,:,:);
% 
% 
% if isfield(SeqParamsUpdated, 'KspaceOrderNavKyKzt')
%     img_nav = ifft3call(NavFully);
%     %     img_nav = img_nav(1+end/4:3*end/4,:,:,:,:);
% end
% 
% 
% for t = 1:num_echoes
%     imagesc3d2(rsos(img(:,:,:,:,t),4), s(img)/2, t, [-90,180,180], [0.,1e-3]), setGcf(.5)
% end

R_factor = zeros(num_echoes,1);

for t = 1:num_echoes
    tmp = Data(:,:,:,1,t);
    R_factor(t) = 1/mean(tmp(:)~=0);
end

disp(R_factor)

% save('data/DataHalfwaythrough.mat', 'img' , 'img_nav', 'B0_smooth', 'Nav_n')
% save('data/DataHalfwaythrough_B0corr.mat', 'img' , 'img_nav', 'B0_smooth')
% save('data/DataHalfwaythrough.mat',  'm2d','-append')
% save('data/DataHalfwaythrough_B0corr.mat',  'm2d','-append')


save([Input.data_path,  Input.DataOutputFilename], 'Data' ,'m2d','-v7.3')
clear Data
clear m2d

if isfield(SeqParamsUpdated, 'KspaceOrderNavKyKzt')
    save([Input.data_path,  Input.DataOutputFilename] , 'NavFully', 'Nav_n','-append','-v7.3')
clear NavFully
clear Nav_n

end;

if Options.Apply_B0_corr
    save([Input.data_path,  Input.DataOutputFilename] , 'B0_smooth','-append','-v7.3')
%     save([Input.data_path,  Input.DataOutputFilename] ,  'b','-append','-v7.3')

end;


%%





