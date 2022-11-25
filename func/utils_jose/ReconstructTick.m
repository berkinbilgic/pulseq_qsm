function img_res = ReconstructTick(Input, Options)
%
% Input.data_path,
% Input.ReceiveFilename
% Input.ReceiveName ='??eceive'
% Input.DataFilename
% Input.DataName ='Data'
% Input.NavFilename
% Input.NavName ='NavFully'
% Input.NavLabeln ='Nav_n'
%
%
% Options.lsqr_iter = 200;
% Options.lsqr_tol = 1e-3;
% Options.lambda = 1e-4;
% Options.useNavdata = 1;
% Options.slices =  ???
% keyboard
receive = load([Input.data_path,Input.ReceiveFilename],Input.ReceiveName);
receive = single(receive.(Input.ReceiveName ));
Data = load([Input.data_path,Input.DataFilename],Input.DataName );
Data = single(Data.(Input.DataName));
% load([Input.data_path, Input.SeqFilename]);

Dims = size(Data);
num_chan = Dims(4);
num_echoes = Dims(5);


if Options.useNavdata
    NavFully = load([Input.data_path,Input.NavFilename],Input.NavName );
    NavFully = NavFully.(Input.NavName);
    Nav_n = load([Input.data_path,Input.NavFilename],Input.NavLabeln );
    Nav_n = Nav_n.(Input.NavLabeln);
    %     load([Input.data_path, Input.NavFilename], 'NavFully', 'Nav_n')
    mask_data = (Data(round(end/2),:,:,1,:)) ~= 0;
    mask_nav = (NavFully(round(end/2),:,:,1,:)) ~= 0;
    Data = (Data + NavFully)./(mask_nav + mask_data + eps);
%     Data (Data==0) = single(NavFully (Data==0));
end



% lsqr_iter = 200;
% lsqr_tol = 1e-3;
% keyboard
for k = 1:num_chan
    for l = 1:num_echoes
        Data(:,:,:,k,l) = fftshift(ifft(ifftshift(Data(:,:,:,k,l), 1), [], 1), 1) * sqrt(Dims(1)); % fourier transform on first dim
    end;
end;
% mask =  squeeze(Data(round(end/2),:,:,1,:)) ;
m2d = squeeze(Data(round(end/2),:,:,:,:)) ~= 0;

img_res = zeross([length(Options.slices),Dims(2:3),num_echoes]);
% mosaic(sq(m2d(:,:,1,:)),2,4,56,'',[0,1]),colormap parula
% mosaic(mean(m2d(:,:,1,:),4),1,1,57,'',[0,1]),colormap parula
% figure(56), imab(sq(m2d(:,:,1,:)),[0,1]),colormap parula
% figure(57), imab(mean(m2d(:,:,1,:),4),[0,1]),colormap parula; colorbar
%



tic
parfor slc_select = Options.slices; %1:s(img,1)
    disp(num2str(slc_select))
    echobyecho = 0 ;
    %     if echobyecho
    %         sens = sq(receive(slc_select,:,:,:));
    %         kspace_slc = (sq(Data(slc_select,:,:,:,:)));
    %
    %
    %         param = [];
    %         param.N = Dims(2:3);
    %         param.num_chan = num_chan;
    %         param.lambda = Options.lambda;        % L2 reg
    %
    %         param.sens = sens;
    %
    %         Res = zeross([param.N, num_echoes]);
    %
    %         for t = 1:num_echoes
    %             kspace_coils = kspace_slc(:,:,:,t);
    %             param.m2d = m2d(:,:,:,t);
    %
    %             res = lsqr(@apply_sense_tikc, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);
    %
    %             Res(:,:,t) = reshape(res, param.N);
    %         end
    %     else
    sens = squeeze(receive(slc_select,:,:,:));
    kspace_slc = (squeeze(Data(slc_select,:,:,:,:)));
    kspace_coils = permute(kspace_slc(:,:,:,:), [1 2 4 3]);
    param = [];
    param.m2d = permute(m2d(:,:,:,:) , [1 2 4 3]);
    param.N = [Dims(2), Dims(3), num_echoes];
    param.sens = permute(sens, [1 2 4 3]);
    param.regul = 0; % uses standard Tikono L2 recon if 0 and inverted recon if weight 1
    param.num_chan = num_chan;
    
    param.lambda = Options.lambda;        % L2 reg
    res = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);
    
    Res = reshape(res, param.N);
    
    
    %     end
    
    img_res(slc_select,:,:,:) = Res;
    %     keyboard
end
toc
img_res = img_res(Options.slices,:,:,:);
% imab(squeeze(img_res(slc_select,:,:,:)))


