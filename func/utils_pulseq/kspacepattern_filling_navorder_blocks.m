function [KspaceNavMaskKyKzt , navcount] = kspacepattern_filling_navorder_blocks(dims,pat)
% dims are the dimensions of the of the image to be acquired x,y,z,t
% pat is a kernal sampling pattern in ky,kz as values 0.5 are the position
% associated with shifts
% the echo times are used to speed up the number of navigators sampled
% during the each train


dimspat = [1 size(pat)];

% ensure the imput dims have 4 dimensions
if length(size(dims))==3
    dims(end+1) = 1;
end

% initializing output variables
KspaceNavMaskKyKzt = zeros(dims(2:4));
KspaceNavMaskKyKzt_sampled = zeros(dims(2:4));
navcount = zeros([prod(dims(2:4)),1]);

[shifty,shiftz] = ind2sub(size(pat),find(pat==0.5));
shiftz(end+1) = 1; shifty(end+1) = 1; % add the default entry

R = (length(shiftz));

% for k=1:R
% subplot(floor(sqrt(R)),ceil(R/floor(sqrt(R))),k)
% imab(circshift(pat,[shifty((k)),shiftz((k))])) ; %the last was repeated
% % imab(circshift(pat,[shiftz((k)),shifty((k))]))
% end



ntarget_smallpat = length(find(pat==1));

rep3 = dims(3)/dimspat(3);
rep2 = dims(2)/dimspat(2);


% this shif ensures that a full sampling is obtained of the k-space center

counter = 0;


ntarget = length(find(pat==1))*rep2*rep3;
% shift = @ (totalEchos,echon,run) +mod(totalEchos-echon+1,totalEchos) - run +1;
% shiftsize_last =0;
% I will now create a mask for the full navigator region
for nnav=1:R % this are the number of shifts needed to make the k-space fully sampled
    patmask = circshift(pat,[shifty((nnav)),shiftz((nnav))]); % makes sure we go
    patmask (patmask~=1 )=0;
    
    
     n = 0;
    BlocksOrder = zeros (rep2*rep3,1);
    for k=1:rep2
        for l=1:rep3
            n = n+1;
            if mod(k,2) ==1
                BlocksOrder( (k-1) * rep3 + l) = n;
            else
                BlocksOrder( (k) * rep3 - (l-1)) = n;
            end;
        end
    end
    Nblock = makemosaic(ones([size(patmask),rep2*rep3])   .* reshape(BlocksOrder,[1,1,rep2*rep3]),rep3);

%     keyboard
    
    patmask = repmat(patmask,rep2,rep3);
    
    for k = 1:rep2*rep3
    [posy([1:ntarget_smallpat] + ((k-1)*ntarget_smallpat)),posz([1:ntarget_smallpat] + ((k-1)*ntarget_smallpat))] = ...
                        ind2sub(size(patmask),find(and(patmask==1,Nblock==k)));
    end;
    %         posy = circshift(posy,shiftsize);
    %         posz = circshift(posz,shiftsize);
    
    shotnumber=zeros(dims(4),ntarget);
    % matrix nruns of navigator vs ntargets  to be filled with
    % Excitation number
    echonumber=zeros(dims(4),ntarget);
    % matrix nruns of navigator vs ntargets  to be filled with
    % echo number. Rule of this matrix there can't be to equal echo
    % numbers in one column
    %     keyboard
    echomatrix = reshape(repmat((1:dims(4)),[1 ntarget]),[ntarget,dims(4)])';
    shotmatrix = reshape(reshape(repmat((1:ntarget),[dims(4),1]),[dims(4)*ntarget 1]),[ntarget,dims(4)])';
    
    % because the number of measurements to obtain a navigator might not be
    % dividible by the number of echo times, there could be very large
    % k-space jumps between echos. By measuring successive navigators in
    % separate direction this problem is reduced
    %
    echomatrix(2:2:end,:) = flipdim(echomatrix(2:2:end,:),1);
    shotmatrix(2:2:end,:) = flipdim(shotmatrix(2:2:end,:),1);
    shotmatrix = shotmatrix + (nnav-1) * ntarget;
    for k = 1:dims(4)
        for l = 1:ntarget
            current_echo = echomatrix(k,l);
            current_shot = shotmatrix(k,l);
            success = 0;
            %shotnumber=zeros(dims(4),ntarget);
            % matrix nruns of navigator vs ntargets  to be filled with
            % Excitation number
            %echonumber=zeros(dims(4),ntarget);
            % matrix nruns of navigator vs ntargets  to be filled with
            % echo number. Rule of this matrix there can't be to equal echo
            % numbers in one column
            availpos = find(shotnumber(k,:)==0);
            pos = 0;
            while success ==0
                pos = pos + 1 ;
                
                if isempty(find(echonumber(:,availpos(pos)) == current_echo))
                    success = 1;
                    echonumber(k,availpos(pos)) = current_echo;
                    shotnumber(k,availpos(pos)) = current_shot;
                    KspaceNavMaskKyKzt( posy(availpos(pos)),posz(availpos(pos)), current_echo) = current_shot; % this is the excitation number
                    KspaceNavMaskKyKzt_sampled(posy(availpos(pos)),posz(availpos(pos)), current_echo) = ...
                        KspaceNavMaskKyKzt_sampled( posy(availpos(pos)),posz(availpos(pos)), current_echo) +1; % this is the excitation number
                    navcount ((current_shot-1)*dims(4) + current_echo) = (nnav -1)*dims(4)+ k;
                    
                end
            end;
        end
    end
    
    
end;

function imall=makemosaic(im,MaxN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function imall=makemosaic(im,MaxN);
% Make a mosaic image for display the image with "show.m"
% i.e., the 3D image [im] transforms to a mosaic 2D image [imall]
% If [im] is 4D, [im(:,:,:,1)] will be used
% NOTE : First, singleton dimensions will be removed; 64x64x1x20 -> 3D
% Input :
%   [im] : 3D or 4D image
%   [MaxN](option): The number of colons, Default is 5
% Output :
%   [imall]: mosaic 2D image
% Usages,
% imall=makemosaic(im,MaxN);
% imall=makemosaic(im);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples,
% let [im] is a matrix of 64x64x20
% imall=makemosaic(im,10);
% [imall] is a 2x10 image of 64x64, size(imall)= 128 x 640
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) Jaemin Shin
% jaemins@gatech.edu
% 01/30/06
% updated 
% 02/28/06 : bug fixed

if exist('MaxN','var') == 0
    MaxN = 5;
end
im = squeeze(im);
dim = size(im);
if length(dim) < 2;
    error('Input is 1D or 2D signal')
elseif length(dim) ==4
    im = squeeze(im(:,:,:,1));
    disp('4D : TimePoint 1 was used')
elseif length(dim) > 4
    error('5D or Higher dimension does not support')
end
Nrow = ceil(dim(3)/MaxN);
Rcol = mod(MaxN - mod(dim(3),MaxN),MaxN);

if dim(3) <= MaxN
    imall = reshape(im,[dim(1) dim(2)*dim(3)]);
    imall = [imall,zeros(dim(1),dim(2)*Rcol)];
else
    imall = reshape(im(:,:,1:MaxN),[dim(1) dim(2)*MaxN]);
    for ii=2:Nrow-1 % bug fixed
        temp = reshape(im(:,:,(ii-1)*MaxN+1:ii*MaxN),[dim(1) dim(2)*MaxN]);
        imall = cat(1,imall,temp);
    end
    temp = reshape(im(:,:,(Nrow-1)*MaxN+1:end),[dim(1) dim(2)*(MaxN-Rcol)]);
    temp = [temp,zeros(dim(1),dim(2)*Rcol)];
    imall = cat(1,imall,temp);
end
