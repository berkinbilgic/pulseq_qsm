function [ res, tflag ] = apply_sense_tikc_JPM( in, params, tflag )
% this function can be used in lsqr to do a sense recosntruction
% params.lambda regularization weigth
% params.regul
%   sense reconstruction of 3D matrices where the either a standard tikonov
%   params.regul == 0
%   or using a tv regularization over the 3rd dimension params.regul == 1
%
%   params.weight are a weightint term used prior to applying the TV
%   regularization



if strcmp(tflag,'transp')
%     display(['transp'])
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels

    b = in(1+params.num_chan*prod(params.N):end);

    kspace_coils = reshape(in(1:(params.num_chan*prod(params.N))), [params.N, params.num_chan]);
    %     display(['kspace_coils', num2str( size(kspace_coils))] )

    img_coils = ifft2call( kspace_coils .* params.m2d );
    %     display(['img_coils',  num2str(size(img_coils))] )
    %     display(['img_coils .* conj(params.sens)',  num2str(size(img_coils .* conj(params.sens)))] )

    %     Res = sum(img_coils .* conj(params.sens), length(size(params.sens))  );
    Res = sum(img_coils .* conj(params.sens), 4  );

    if params.regul == 1 % weighted gradient along z
%         b = reshape(in(1+params.num_chan*prod(params.N):end),[params.N params.num_chan]);
        b = reshape(in(1+params.num_chan*prod(params.N):end),[params.N ]);
        b = conj(params.weight) .* gradz( b ,'transp',[],-1);
    end


    res = Res(:) + sqrt(params.lambda) * b(:);

else
    % display(['notransp'])

    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    in = reshape(in, params.N);
    img_coils = repmat(in, [1,1,1,params.num_chan]) .* params.sens;
    %     display(['img_coils',  num2str(size(img_coils))] )
    %     display(['fft2call(img_coils)',  num2str(size(fft2call(img_coils)))] )
    kspace_coils = fft2call(img_coils) .* params.m2d;
    %     display(['kspace_coils',  num2str(size(kspace_coils))] )

    if params.regul == 1 % weighted gradient along z
        in = gradz(params.weight .* in ,'notransp',[],-1);
    end
    res = cat(1, kspace_coils(:), sqrt(params.lambda)*in(:));

end

end
