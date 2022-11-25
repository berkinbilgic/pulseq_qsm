function [ img_res, img_sum] = ReconstructTick_nav(Input, Options)
%
% Input.data_path,
% Input.ReceiveFilename
% Input.ReceiveName ='ŕeceive'
% Input.NavFilename
% Input.NavName ='NavFully'
% Input.NavLabeln ='Nav_n'
%
%
% Options.lsqr_iter = 200;
% Options.lsqr_tol = 1e-3;
% Options.lambda = 1e-4;
% Options.nav_numb =  ???
% keyboard
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


%keyboard

% lsqr_iter = 200;
% lsqr_tol = 1e-3;
Dims = size(Data);

% if isfield(Options,'Nx')
% Dims (1) = Options,'Nx';
% end
Data = fftshift(ifft(ifftshift(Data, 1), [], 1), 1) * sqrt(Dims(1)); % fourier transform on first dim
% Data = fftshift(ifft(ifftshift(crop(Data,[ Dims(2:end)]), 1), [], 1), 1) * sqrt(Dims(1)); % fourier transform on first dim
% mask =  squeeze(Data(round(end/2),:,:,1,:)) ;
m2d = squeeze(Data(round(end/2),:,:,:,:)) ~= 0;

num_chan = Dims(4);
num_echoes = Dims(5)
img_res = zeross([Dims(1:3),num_echoes]);
% mosaic(sq(m2d(:,:,1,:)),2,4,56,'',[0,1]),colormap parula
% mosaic(mean(m2d(:,:,1,:),4),1,1,57,'',[0,1]),colormap parula
% figure(56), imab(sq(m2d(:,:,1,:)),[0,1]),colormap parula
% figure(57), imab(mean(m2d(:,:,1,:),4),[0,1]),colormap parula; colorbar
%


% delete(gcp('nocreate'))
% c = parcluster('local');    % build the 'local' cluster object

% total_cores = c.NumWorkers;
% parpool(ceil(total_cores * .5))

tic
for slc_select = 1:Dims(1)
    disp(num2str(slc_select))

    sens = sq(receive(slc_select,:,:,:));
    kspace_slc = (sq(Data(slc_select,:,:,:,:)));


    param = [];
    param.N = Dims(2:3);
    param.num_chan = num_chan;
    param.lambda = Options.lambda;        % L2 reg

    param.sens = sens;

    Res = zeross([param.N, num_echoes]);

    for t = 1:num_echoes
        kspace_coils = kspace_slc(:,:,:,t);
        param.m2d = m2d(:,:,:,t);

        res = lsqr(@apply_sense_tikc, cat(1, kspace_coils(:), zeros(prod(param.N),1)), Options.lsqr_tol, Options.lsqr_iter, [], [], [], param);

        Res(:,:,t) = reshape(res, param.N);
    end
    img_res(slc_select,:,:,:) = Res;

end

% if length(Options.slices)==1
%     img_res = Res;
% end;
% 
% imab(squeeze(img_res(slc_select,:,:,:)))
img_sum = sum(abs(img_res),4);
% delete(gcp('nocreate'))


% for t = 1:num_echoes
%     imagesc3d2(img_res(:,:,:,t), s(img_res)/2, t, [-90,180,-180], [0.,3e-3]), setGcf(.5)
% end
