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
% Options.nav_numb =  ???
% currently this is not actually being used
% Options.regul=  ???% uses standard Tikono L2 recon if 0 and time inverted weight if weight 1
% Options.weight=  ???% time inverted estimation



receive = load([Input.data_path,Input.ReceiveFilename],Input.ReceiveName);
receive = receive.(Input.ReceiveName );

% load([Input.data_path, Input.SeqFilename]);


NavFully = load([Input.data_path,Input.NavFilename],Input.NavName );
NavFully = NavFully.(Input.NavName);
Nav_n = load([Input.data_path,Input.NavFilename],Input.NavLabeln );
Nav_n = Nav_n.(Input.NavLabeln);
%     load([Input.data_path, Input.NavFilename], 'NavFully', 'Nav_n')



mask = 0*Nav_n ;
for mc = Options.nav_numb
    mask = double(Nav_n == mc) +mask;
end

Data = NavFully .* double(mask);


Dims = size(Data);
num_chan = Dims(4);
num_echoes = Dims(5)

for k = 1:num_chan
    for l = 1:num_echoes
        Data(:,:,:,k,l) = fftshift(ifft(ifftshift(Data(:,:,:,k,l), 1), [], 1), 1) * sqrt(Dims(1)); % fourier transform on first dim
    end;
end;

m2d = squeeze(Data(round(end/2),:,:,:,:)) ~= 0;

img_res = zeross([Dims(1:3),num_echoes]);


% keyboard
% tic
for slc_select = 1:Dims(1)
    disp(num2str(slc_select))
    
    sens = sq(receive(slc_select,:,:,:));
    
    
    
    kspace_slc = (sq(Data(slc_select,:,:,:,:)));
    kspace_coils = permute(kspace_slc(:,:,:,:), [1 2 4 3]);
    param.m2d = permute(m2d(:,:,:,:) , [1 2 4 3]);
    param.N = [Dims(2), Dims(3), num_echoes];
    param.sens = permute(sens, [1 2 4 3]);
    param.regul = Options.regul; % uses standard Tikonov L2 recon if 0 and inverted recon if weight 1
    param.num_chan = num_chan;
    
    param.lambda = Options.lambda;        % L2 reg
    %     if isfield(Options,'weights')
    param.weight = squeeze(Options.weight(slc_select,:,:,:));
    res = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(Dims(2)* Dims(3)*num_echoes,1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);
    
    %     else
    %
    %     param.weight=flipdim((reshape(res2, param.N)),3);
    %     res2 = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], res2, param);
    
    %
    %     end;
    img_res(slc_select,:,:,:) = reshape(res, param.N);
    
    % just for debugging purposes
    % param.lambda  = 1
    % res = lsqr(@apply_sense_tikc_JPM, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);
    % imab(reshape(res, param.N))
    
end

img_sum = sum(abs(img_res),4);

