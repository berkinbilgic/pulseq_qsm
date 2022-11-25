function MakeGif (filename,images_cellarray,title_cellarray,range,varargin)

if nargin < 4 || isempty(range)
    range = prctile(images_cellarray{1}(:),[2 98]);
end
if nargin == 5
    saveseparatily = varargin{1};
else
    
    saveseparatily =0;
    
end;

figure(1)

set(gcf,'Color',[1 1 1]);

for n = 1:length(images_cellarray)
    imagesc(images_cellarray{n},range)
    colormap(gray)
    axis equal
    axis off
    axis tight
    colorbar
    if nargin < 3 || isempty(title_cellarray)
        title([ 'frame ',num2str(n)])
    else
        title(title_cellarray{n}, 'fontSize',12);
    end
    drawnow
    frame = getframe(1);
    im (:,:,:,n)= frame2im(frame);
    
    %       [imind,cm] = rgb2ind(im,256);
end

if n == 1;
    imwrite(im(:,:,1,:),filename,'gif', 'Loopcount',inf,'DelayTime',1);
else
    %           imwrite(im,filename,'gif','WriteMode','append','DelayTime',1);
    %           imwrite(im,filename,'gif', 'Loopcount',3,'DelayTime',1);
    %           imwrite(im,'testloopinf.gif','gif', 'Loopcount',inf,'DelayTime',1);
    %           imwrite(im(:,:,1,:),'testloopgrey.gif','gif', 'Loopcount',inf,'DelayTime',1);
    if saveseparatily == 1
        for count = 1: n
                    display('I amhere')

            imwrite(im(:,:,1,count),[filename,'_frame_',num2str(count), '.jpg'],'jpeg','Quality',90);
            
        end
    else
    imwrite(im(:,:,1,:),filename,'gif', 'Loopcount',inf,'DelayTime',1);
        
    end
end
