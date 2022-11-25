function [parameters, NRMSE, l2loss_DC]  = fitting_S0_R2_chi_NativeRes_tmp_with_BackgroundFieldUpsampling(ss,te,mask,parameters,FitOptions);
% ss 4d multiecho data to be fitted
% te
% D = B0 * gyro * create_dipole_kernel(B0_dir , voxelSize, matrixSize, 1);
% mask original is the brain mask where the fitting will be performed
% parameters use a structure with dlarrays of the parameters to be fitted...
% this
% FitOptions = struct('numEpochs', 400, 'fitonlymask', 1, 'initialLearnRate', 0.1,'decayRate', 0.0001);
% gpu = gpuDevice();

if isempty(FitOptions)
    % Specify Training Options
    numEpochs = 400;
    % executionEnvironment = "gpu";
    initialLearnRate = 0.1;
    decayRate = 0.0001;
    fitonlymask =true;
    %     SpatialUnwrapEpoch = 0;
    %     FitOptions.sTVregul = 0; %default is no STV regul
    %     FitOptions.sTVlambda = 1;  % regularization weigth
    %     FitOptions.sTVweight = [];  % regularization weigth
    %     draw = false;
    %     FitOptions.magn_phase_chi_fit = ['110'] ;% defined what is fitted
    %     FitOptions.NormalizeL2loss = 1;
    %     FitOptions.Field_5D = 0;

else
    if isfield(FitOptions,'numEpochs')
        numEpochs = FitOptions.numEpochs;
    else
        numEpochs = 400;
    end
    %     if isfield(FitOptions,'SpatialUnwrapEpoch')
    %         SpatialUnwrapEpoch = FitOptions.SpatialUnwrapEpoch;
    %     else
    %         SpatialUnwrapEpoch = numEpochs / 2;
    %     end
    if isfield(FitOptions,'initialLearnRate')
        initialLearnRate = FitOptions.initialLearnRate;
    else
        initialLearnRate = 0.1;
    end
    if isfield(FitOptions,'decayRate')
        decayRate = FitOptions.decayRate;
    else
        decayRate = 0.0001;
    end
    if isfield(FitOptions,'fitonlymask')
        fitonlymask = FitOptions.fitonlymask;
    else
        fitonlymask = 1;
    end
    if isfield(FitOptions,'sense')
        FitOptions.sense = FitOptions.sense;
    else
        FitOptions.sense = 1;
    end
    %     if isfield(FitOptions,'sTVlambda')
    %         FitOptions.sTVlambda = FitOptions.sTVlambda;
    %     else
    %         FitOptions.sTVlambda = 0;
    %     end
    %     if ~isfield(FitOptions,'sTVweight')
    %         FitOptions.sTVweight = [];
    %     end
    %
    %     if ~isfield(FitOptions,'D')
    %         FitOptions.D = [];
    %     end
    if isfield(FitOptions,'draw')
        draw = FitOptions.draw;
    else
        draw = true;
    end
    %     if isfield(FitOptions,'magn_phase_chi_fit')
    %         FitOptions.magn_phase_chi_fit  = FitOptions.magn_phase_chi_fit ;
    %         %         'should be here'
    %     else
    %         'ups'
    %         FitOptions.magn_phase_chi_fit  = '110';
    %     end
    %     if ~isfield(FitOptions,'NormalizeL2loss')
    %         FitOptions.NormalizeL2loss = 0;
    %     end
    %
    %     % this could also be called a background field
    %     if isfield(FitOptions,'Field_5D')
    %         FitOptions.Field_5D;
    %     else
    %         FitOptions.Field_5D = 0;
    %     end
    %

end
% prepare input data
[sx, sy, sz, st] = size(ss);
if isempty(mask)
    mask = ones([sx, sy, sz]);
end;

scaleFactor = sqrt(sum(abs(ss).^2, 4))/size(FitOptions.Field_5D,5);


% keyboard
if fitonlymask
    scales = struct('s0', scaleFactor.*mask, 'r2s', 50,'freq', 20,  'phi0', 10); % this parameters are related to the amplitude of the metrics

else
    scales = struct('s0', scaleFactor, 'r2s', 50, 'freq', 20, 'phi0', 10); % this parameters are related to the amplitude of the metrics
end

dTE = mean(diff(te));
te = reshape (te,[1 1 1 st]);

% initialize model parameters

% parameters = initialise_model_r2s(ss(:,:,:,1));

if isempty(parameters)
    parameters = struct;

    parameters.s0 = dlarray(rand(sx, sy, sz, 'single'));
    parameters.r2s = dlarray(rand(sx, sy, sz, 'single')) * 10;
    parameters.phi0 = dlarray(rand(sx, sy, sz, 'single'))* 0;
    parameters.freq = dlarray(rand(sx, sy, sz, 'single')) * 0;

else
    parameters.s0=dlarray(parameters.s0./scales.s0);
    parameters.s0(isnan(parameters.s0))=0;
    parameters.s0(isinf(parameters.s0))=0;

    parameters.r2s=dlarray(parameters.r2s./scales.r2s);
    parameters.r2s(isnan(parameters.r2s))=0;
    parameters.r2s(isinf(parameters.r2s))=0;

    parameters.phi0=dlarray(parameters.phi0./scales.phi0);
    parameters.phi0(isnan(parameters.phi0))=0;
    parameters.phi0(isinf(parameters.phi0))=0;

    parameters.freq=dlarray(parameters.freq./scales.freq);
    parameters.freq(isnan(parameters.freq))=0;
    parameters.freq(isinf(parameters.freq))=0;

end
display (' creating GPU arrays')
parameters.s0=gpuArray(parameters.s0);
parameters.r2s=gpuArray(parameters.r2s);
parameters.phi0=gpuArray(parameters.phi0);
parameters.freq=gpuArray(parameters.freq);


% keyboard
accfun = dlaccelerate(@modelGradients_r2s);
% keyboard
% prepare data for gpu
scales.ind = gpuArray(find(((ones(size(ss)).*mask))~=0));
mask = gpuArray(dlarray(mask));


ss_real = real(ss);
ss_imag = imag(ss);
ss_real = gpuArray(dlarray(ss_real));
ss_imag = gpuArray(dlarray(ss_imag));
if fitonlymask
    display(' Data outside Mask = 0')

    ss_real = ss_real .*mask;
    ss_imag = ss_imag .*mask;
end
if FitOptions.NormalizeL2loss
    FitOptions.NormalizeL2loss = (sum(abs(ss(scales.ind )).^2))/numel(scales.ind)/2;
else
    FitOptions.NormalizeL2loss = 1;
end
clear ss

te = gpuArray(dlarray(te));
% mask = gpuArray(dlarray(repmat(maskOriginal, [1 1 1 length(te)])));

% scales.ind = gpuArray(find(extractdata(gather(ss_real))~=0));
scales.ind = gpuArray(find(extractdata(gather(ones(size(ss_real)).*mask))~=0));
scales.s0 = gpuArray(dlarray(scales.s0));


% if ~isempty(FitOptions.sTVweight)
%     % aim is to discontinue the scales
%     scales.sTVweight=gpuArray(dlarray(FitOptions.sTVweight));
%     FitOptions.sTVweight=gpuArray(dlarray(FitOptions.sTVweight));
% end;

if ~isempty(FitOptions.sense)
    % aim is to discontinue the scales
    FitOptions.sense=gpuArray(dlarray(FitOptions.sense));
end;


averageGrad = [];
averageSqGrad = [];
% keyboard
figure;
a = gcf;
C = colororder;
lineLoss = animatedline('Color', C(2, :));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

start = tic;

iteration = 0;
scounter = 1;

generate_plots(parameters, mask,scales);

for epoch = 1:numEpochs
    iteration = iteration + 1;
    %      keyboard
    % Evaluate the model gradients and loss using dlfeval and the
    % modelGradients function.

    [gradients, loss, l2loss_DC] = dlfeval(accfun, parameters, ss_real, ss_imag, te, FitOptions, scales , freq);

    % Update learning rate.
    learningRate = initialLearnRate / (1 + decayRate * iteration);

    % Update the network parameters using the adamupdate function.
    [parameters, averageGrad, averageSqGrad] = adamupdate(parameters, gradients, averageGrad, ...
        averageSqGrad, iteration, learningRate);


    %         make sure no nan appear
    if strcmp(FitOptions.magn_phase_chi_fit(1),'1')

        parameters.s0( isnan(parameters.s0) ) = 0;

        parameters.r2s( isnan(parameters.r2s) ) = 0;

    end
    if strcmp(FitOptions.magn_phase_chi_fit(2),'1')

        parameters.phi0( isnan(parameters.phi0) ) = 0;

        if strcmp(FitOptions.magn_phase_chi_fit(3),'1')

            parameters.chi( isnan(parameters.chi) ) = 0;

        else
            parameters.freq( isnan(parameters.freq) ) = 0;

        end

    end;

    % Plot training progress.
    loss = double(gather(extractdata(loss)));
    addpoints(lineLoss, iteration, loss);

    if and (draw ,mod(epoch, 10) == 0)
        figure(a)
        Dur = duration(0, 0, toc(start), 'Format', 'hh:mm:ss');
        title("Epoch: " + epoch + ", Elapsed: " + string(Dur) + ", Loss: " + loss)
        drawnow
    end
    if mod(epoch, 100) == 0
        % print intermediate loss
        %         [intm_real, intm_imag] = model_r2s_chi(parameters, te, FitOptions.D, scales,[]);


        [R_real, R_imag] = model_complex(parameters, te, FitOptions, scales,freq);
        d_real = gather(extractdata(abs(ss_real - R_real))).^2;
        d_imag = gather(extractdata(abs(ss_imag - R_imag))).^2;
        d_total = d_real + d_imag;
        NRSE = sqrt(sum((d_total .* mask), 'all')) / sqrt(sum((ss_real.^2 + ss_imag.^2) .* mask, 'all'));

        % keyboard
        %         NRSE = sqrt(sum((d_total .* mask), 'all')) / sqrt(sum((ss_real.^2 + ss_imag.^2) .* mask, 'all'));
        fprintf('#Epoch: %i, NRSE: %f \n', epoch, NRSE);
    end


end
gpu = gpuDevice();
display(['gpu.AvailableMemory ', gpu.AvailableMemory]);

showReults = 1;
if showReults
    [R_real, R_imag] = model_complex(parameters, te, FitOptions, scales,freq);
    d_real = gather(extractdata(abs(ss_real - R_real))).^2;
    d_imag = gather(extractdata(abs(ss_imag - R_imag))).^2;
    d_total = d_real + d_imag;

    NRMSE_show = sqrt(sum(d_total, 4))./scales.s0;
    NRMSE = NRMSE_show ;
    NRMSE =  gather(extractdata(NRMSE));
    NRMSE_show(NRMSE <= 0.05 & NRMSE ~= 0) = 0.5;
    NRMSE_show(NRMSE > 0.05) = 1;

    figure(45);
    subplot(2, 1, 1);
    Orthoview(NRMSE .* mask, [], [0, 1], [], jet);
    title('NRMSE explanation');
    colorbar;

    subplot(2, 1, 2);
    Orthoview(NRMSE_show .* mask, [], [0, 1], [], cat(1, [0 0 0], jet));
    title('NRMSE explanation (bellow 0.05)');

    generate_plots(parameters, mask, scales);

end


parameters.phi0 = angle(exp(i*(gather(extractdata(parameters.phi0)) * scales.phi0 + angle(gather(extractdata(parameters.s0))))));
parameters.freq = gather(extractdata(parameters.freq .* scales.freq));
parameters.r2s = gather(extractdata(parameters.r2s .* scales.r2s));
parameters.s0 = abs(gather(extractdata(parameters.s0 .* scales.s0)));
l2loss_DC = gather(extractdata(l2loss_DC));
end

function [gradients, loss, l2loss_DC] = modelGradients_r2s(parameters, ss_real, ss_imag, te, FitOptions, scales,freq);

% Make predictions with the initial conditions.
% keyboard
[R_real, R_imag] = model_complex(parameters, te, FitOptions, scales,freq);

loss = l2loss([R_real(scales.ind),R_imag(scales.ind)],...
    [ ss_real(scales.ind),ss_imag(scales.ind)],...
    'DataFormat','CB','NormalizationFactor','all-elements')/FitOptions.NormalizeL2loss;

l2loss_DC = loss;

if isfield(FitOptions,'sTVweight')
    if size(FitOptions.sTVweight,4)==3

        sTVhandle = str2func('sTVmat_fast');
    else
        if size(FitOptions.sTVweight,4)==6
            sTVhandle = str2func('sTVmat');
        else
            if isempty(FitOptions.sTVweight)
                sTVhandle = str2func('sTVmat_fast');
            end
        end
    end;
end;

if or (FitOptions.sTVregul == 2,FitOptions.sTVregul == 1)

    % if chi is being fitted the regularization is on chi, else it is on
    % r2*
    if strcmp(FitOptions.magn_phase_chi_fit(3),'1')
        if isempty(FitOptions.sTVweight)
            sTVregul = sTVhandle(parameters.chi,parameters.r2s.*parameters.s0 , scales.s0~=0);

        else
            sTVregul = sTVhandle(parameters.chi,FitOptions.sTVweight , scales.s0~=0);
        end
    else
        if isempty(FitOptions.sTVweight)
            sTVregul = sTVhandle(parameters.r2s, parameters.s0 , scales.s0~=0);
        else
            sTVregul = sTVhandle(parameters.r2s, FitOptions.sTVweight , scales.s0~=0);
        end
    end

    if FitOptions.sTVregul == 2

        loss = loss + l2loss(FitOptions.sTVlambda * sTVregul(:), ...
            0*sTVregul(:),...
            'DataFormat','CB','NormalizationFactor','all-elements');
    end
    if FitOptions.sTVregul == 1

        loss = loss + l1loss(FitOptions.sTVlambda * sTVregul(:), ...
            0*sTVregul(:),...
            'DataFormat','CB','NormalizationFactor','all-elements');
    end
end

% if FitOptions.sTVregul == 1
%
%     % if chi is being fitted the regularization is on chi, else it is on
%     % r2*
%     if strcmp(FitOptions.magn_phase_chi_fit(3),'1')
%         sTVregul = sTVhandle(parameters.chi, FitOptions.sTVweight , scales.s0~=0);
%     else
%         sTVregul = sTVhandle(parameters.r2s, FitOptions.sTVweight , scales.s0~=0);
%     end
%     loss = loss + l1loss(FitOptions.sTVlambda * sTVregul(:), ...
%         0*sTVregul(:),...
%         'DataFormat','CB','NormalizationFactor','all-elements');
%
% end
% keyboard
% Calculate gradients with respect to the learnable parameters.
gradients = dlgradient(loss, parameters);
%
%
if strcmp(FitOptions.magn_phase_chi_fit(2),'1')
    %     gradients.phi0 = real(gradients.phi0);
    if strcmp(FitOptions.magn_phase_chi_fit(3),'1')

        gradients.chi = real(gradients.chi);
        %         gradients.chi (isnan(gradients.chi)) = 0;
        %          gradients.chi (isinf(gradients.chi)) = 0;
        %     else
        %         gradients.freq = real(gradients.freq);

    end
end
if strcmp(FitOptions.magn_phase_chi_fit(1),'1')
    gradients.s0 = real(gradients.s0);
    gradients.r2s = real(gradients.r2s);
end
%  keyboard
end




function [dlU_Real, dlU_Imag] = model_complex(parameters, te, FitOptions, scales , freq)
% freq is only passed to make sure that we are not recreating frequency variable every
% single time
if nargin < 3
    scales = struct();
end

if ~isfield(scales, 'r2s')
    scales.r2s = 1;
end

if ~isfield(scales, 'chi')
    scales.chi = 1;
end

if ~isfield(scales, 'phi0')
    scales.phi0 = 1;
end
if ~isfield(scales, 's0')
    scales.s0 = 1;
end


dlU_Real = (real(scales.s0 .* parameters.s0 .* exp(-te .* parameters.r2s .* scales.r2s) .* cos((parameters.phi0 .* scales.phi0) + 2 .* pi .* (te .* (parameters.freq * scales.freq  )))));

dlU_Imag = (real(scales.s0 .* parameters.s0 .* exp(-te .* parameters.r2s .* scales.r2s) .* sin((parameters.phi0 .* scales.phi0) + 2 .* pi .* (te .* (parameters.freq * scales.freq )))));


end

function out = fft3s_(data)
out = fft(fft(fft(data,[],1),[],2),[],3);
end

function out = ifft3s_(data)
out = ifft(ifft(ifft(data,[],1),[],2),[],3);
end

function sTVnorm  = sTVmat_fast (input, weigth, mask)
% This function computes the sTVnorm using a precomputed sTVweight and applies it to chi
% the output is a 4D matrix with 3 elements on the 4th dimension (that penalize gradients in directions different than the imaged used to compute the weights)

s = size(input);

grad_input = cat(4,cat(1,input(2:end,:,:)-input(1:end-1,:,:),zeros([1 s(2) s(3)])),...
    cat(2,input(:,2:end,:)-input(:,1:end-1,:),zeros([s(1) 1 s(3)])),...
    cat(3,input(:,:,2:end)-input(:,:,1:end-1),zeros([s(1) s(2) 1])));
if size(weigth,4)==1
    grad_weigth = cat(4,cat(1,weigth(2:end,:,:)-weigth(1:end-1,:,:),zeros([1 s(2) s(3)])),...
        cat(2,weigth(:,2:end,:)-weigth(:,1:end-1,:),zeros([s(1) 1 s(3)])),...
        cat(3,weigth(:,:,2:end)-weigth(:,:,1:end-1),zeros([s(1) s(2) 1])));

    structuralmean = 0.01*sum(sum(sum(sqrt(sum(mask.*grad_weigth.^2,4)),1),2),3)/sum(mask(:));


    norm_s0 = sum(grad_weigth.^2,4)+structuralmean^2;
    weigth = single(grad_weigth./sqrt(norm_s0));


    % both sides will be normalized

    %     structuralmean = 0.01*sum(sum(sum(sqrt(sum(mask.*grad_input.^2,4)),1),2),3)/sum(mask(:));
    %
    %
    %     norm_s0 = sum(grad_input.^2,4)+structuralmean^2;
    %     grad_input = grad_input./sqrt(norm_s0).*mask;

end
sTVnorm = (grad_input - sum(weigth.*grad_input,4).*weigth).*mask;

end

function sTVnorm  = sTVmat (input, weigth, mask)
% This function computes the sTVnorm and has two possible behaviours:
% i) computes the STVweight based on s0 and r2s and applies it to chi
% ii) it uses a precomputed sTVweight and applies it to chi
% the output is a 4D matrix with 3 elements on the 4th dimension (that penalize gradients in directions different than the S0 image)

% [chi_x,chi_y,chi_z] = gradient(parameters.chi);
s = size(input);

grad_input = cat(4,cat(1,input(2:end,:,:)-input(1:end-1,:,:),zeros([1 s(2) s(3)])),...
    cat(2,input(:,2:end,:)-input(:,1:end-1,:),zeros([s(1) 1 s(3)])),...
    cat(3,input(:,:,2:end)-input(:,:,1:end-1),zeros([s(1) s(2) 1])));


sTVnorm = mask .* cat(4, ...
    weigth(:,:,:,1).* grad_input(:,:,:,1)...
    +weigth(:,:,:,2).*grad_input(:,:,:,2)...
    +weigth(:,:,:,3).*grad_input(:,:,:,3),...
    weigth(:,:,:,2).* grad_input(:,:,:,1)...
    +weigth(:,:,:,4).*grad_input(:,:,:,2)...
    +weigth(:,:,:,5).*grad_input(:,:,:,3),...
    weigth(:,:,:,3).* grad_input(:,:,:,1)...
    +weigth(:,:,:,5).*grad_input(:,:,:,2)...
    +weigth(:,:,:,6).*grad_input(:,:,:,3));

end

