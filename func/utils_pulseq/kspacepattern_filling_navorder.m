function [KspaceNavMaskKyKzt , navcount] = kspacepattern_filling_navorder(dims,pat)
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



ntarget = length(find(pat==1));



% this shif ensures that a full sampling is obtained of the k-space center

counter = 0;


ntarget = length(find(pat==1))*dims(2)/dimspat(2)*dims(3)/dimspat(3);
% shift = @ (totalEchos,echon,run) +mod(totalEchos-echon+1,totalEchos) - run +1;
% shiftsize_last =0;
% I will now create a mask for the full navigator region
for nnav=1:R % this are the number of shifts needed to make the k-space fully sampled
    patmask = circshift(pat,[shifty((nnav)),shiftz((nnav))]); % makes sure we go
    patmask (patmask~=1 )=0;
    
    
    patmask = repmat(patmask,dims(2)/dimspat(2),dims(3)/dimspat(3));
    [posy,posz] = ind2sub(size(patmask),find(patmask==1));
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
