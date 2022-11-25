function [ dz] = gradz(I,tflag,varargin)
% tfalg the default (if left empty) is tflag='notransp'
% if tflag='transp' is used it is actually the transpose of the gradient
% that os used data is actually - the gradient and has some particularities
% on the edges
% varargin{1} 
% varargin{2} corresponds to the gradientmethod, the default is 0 which corresponds
% to the mid gradient definition
if or(nargin==1,nargin==2)
    res=[1,1,1];
else
    if isempty(varargin{1})
        res=[1,1,1];
    else
        res=varargin{1};
    end;
end

if nargin>=4
    gradientmethod = varargin{2};
else
    gradientmethod = 0;
end;


gradientmethod;
if strcmp(tflag,'transp')
    
    switch gradientmethod
        case 1
            %         disp('forward one');
            %defined based on div_op by Gilles Puy, the forward gradient
            dz=-cat(3,I( :,:,1,:) , I( :,:,2:end-1,:)-I( :,:,1:end-2,:) , -I( :,:,end-1,:));
        case 0
            %         disp('zero');
            %defined based on my favorite way, the mid gradient, see entry in the labbook "Jose Projects"on the 7/11/2011
            dz=-cat(3,I( :,:,1,:)+0.5*I( :,:,2,:) ,-I( :,:,1,:)+0.5*I( :,:,3,:) , 0.5*(I( :,:,4:end-1,:)-I( :,:,2:end-3,:)) , I( :,:,end,:)-0.5*I( :,:,end-2,:), -I( :,:,end,:)-0.5*I( :,:,end-1,:));
        case -1
            %         disp('backward one');
            %defined based on div_op by Gilles Puy, the backward gradient
            dz=-cat(3,I( :,:,2,:) , I( :,:,3:end,:)-I( :,:,2:end-1,:) , -I( :,:,end,:));
        case 2
            %         keyboard
            disp('forward two.. not completely correct n the edges');
            
            dz=-cat(3,I( :,:,1,:) ,I( :,:,2,:) , I( :,:,3,:), 0.25*(I( :,:,6:end-1,:)-I( :,:,2:end-5,:)) , -I( :,:,end-3,:), -I( :,:,end-2,:), -I( :,:,end-1,:));
            
        otherwise
            disp('divergence not performed');
    end
    
    
    
    
    
    
else
    switch gradientmethod
        case 1
            %defined based on gradient_op by Gilles Puy - it is a forward differential
            dz = cat(3, I(:,:,2:end,:)-I(:,:,1:end-1,:) , zeros(size(I, 1), size(I, 2), 1,size(I, 4)));
        case 0
            %defined based on my favorite way - midpoint differenttial
            dz = cat(3, I(:,:,2,:)-I(:,:,1,:) , 0.5*(I(:,:,3:end,:)-I(:,:,1:end-2,:)) ,I(:,:,end,:)-I(:,:,end-1,:));
        case -1
            %defined based on gradient_op by Gilles Puy - it is a backwards differential
            dz = cat(3, zeros(size(I, 1), size(I, 2), 1,size(I, 4)), I(:,:,2:end,:)-I(:,:,1:end-1,:) );
        case 2
            %          keyboard
            %defined based on my favorite way - but very large midpoint differenttial
            dz = cat(3, I(:,:,2:3,:)-I(:,:,1:2,:) , 0.25*(I(:,:,5:end,:)-I(:,:,1:end-4,:)) ,I(:,:,(end-1):(end),:)-I(:,:,(end-2):(end-1),:));
            
        otherwise
            disp('gradient is not performed');
    end
end;
dz = dz /res(3);

end

