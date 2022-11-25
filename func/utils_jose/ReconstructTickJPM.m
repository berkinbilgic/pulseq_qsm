function [ img_res, img_sum] = ReconstructTickJPM_nav(Input, Options)
%
% Input.data_path,
% Input.ReceiveFilename
% Input.ReceiveName ='??eceive'
% Input.NavFilename
% Input.NavName ='NavFully'
% Input.NavLabeln ='Nav_n'
%
%
% Options.lsqr_iter = 200;
% Options.lsqr_tol = 1e-3;
% Options.lambda = 1e-4;
% Options.regul=  ???% uses standard Tikono L2 recon if 0 and time inverted weight if weight 1
% Options.weight=  ???% tiem inverted estimation
% Options.useNavdata = 1;
% Options.slices =  ???

receive = load([Input.data_path,Input.ReceiveFilename],Input.ReceiveName);
receive = single(receive.(Input.ReceiveName ));
Data = load([Input.data_path,Input.DataFilename],Input.DataName );
Data = single(Data.(Input.DataName));
% load([Input.data_path, Input.SeqFilename]);


if Options.useNavdata
    NavFully = load([Input.data_path,Input.NavFilename],Input.NavName );
    NavFully = NavFully.(Input.NavName);
    Nav_n = load([Input.data_path,Input.NavFilename],Input.NavLabeln );
    Nav_n = Nav_n.(Input.NavLabeln);
    %     load([Input.data_path, Input.NavFilename], 'NavFully', 'Nav_n')
    Data (Data==0) = single(NavFully (Data==0));
end



% lsqr_iter = 200;
% lsqr_tol = 1e-3;
Dims = size(Data);
num_chan = Dims(4);
num_echoes = Dims(5);
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


% keyboard
tic
% parfor slc_select = Options.slices; %1:s(img,1)
parfor slc_select = Options.slices; %1:s(img,1)
    disp(num2str(slc_select))

    sens = sq(receive(slc_select,:,:,:));

    param = [];

    kspace_slc = (sq(Data(slc_select,:,:,:,:)));
    kspace_coils = permute(kspace_slc(:,:,:,:), [1 2 4 3]);
    param.m2d = permute(m2d(:,:,:,:) , [1 2 4 3]);
    param.N = [Dims(2), Dims(3), num_echoes];
    param.sens = permute(sens, [1 2 4 3]);
    param.regul =Options.regul; % uses standard Tikono L2 recon if 0 and inverted recon if weight 1
    param.num_chan = num_chan;
    param.lambda = Options.lambda;        % L2 reg
    if isfield(Options,'weight')
        param.weight = squeeze(Options.weight(slc_select,:,:,:));
        res2 = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);

    else
        res2 = zeros(param.N);
        for k =1 :3
            if max(abs(res2(:)))~=0
                param.weight=flip((reshape(res2, param.N)),3)/(norm(res2(:))/sqrt(numel(find(res2(:)~=0))));

            else
                param.weight =res2;
            end

            res2 = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], res2(:), param);
        end

    end;
    img_res(slc_select,:,:,:) = reshape(res2, param.N);

    % just for debugging purposes
    % param.lambda  = 1
    % res = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);
    % imab(reshape(res, param.N))

end
img_res = img_res(Options.slices,:,:,:);

img_sum = sum(abs(img_res),4);

