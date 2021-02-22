function [] = mkActivityGif(folderpath,ds)

    doColors = [85 205 252; 247 168 184; ...
        255 255 255; 247 168 184; 85 205 252]./255;

    doColors = [1 1 1];
    
    clc

    filePrefix = 'msvideo';
    if nargin < 2 || isempty(ds)
        ds = 5; % downsampling
    end
    
    ap = getFilePaths(folderpath,'.avi');
    for i = length(ap):-1:1
        slind = find(ismember(ap{i},'\/'));
        if ~ismember({ap{i}(slind(end)+1:end-4)},{filePrefix})
            ap(i) = [];
        end
    end
    
    for p = ap'
    
        slind = find(ismember(p{1},'/\'));
        outName = ['Plots/Activity_GIFs/' ...
            p{1}(slind(end-3)+1:slind(end-1)-1) '_' num2str(ds) 'x'];
        
        if exist([outName '_dFF.gif'],'file')==2
            fprintf(['\n\t\t\tAlready created:  ' outName])
            continue
        end
        
        fprintf(['\n\t\t\tCreating:  ' outName])
        
        vidObj{1} = VideoReader(p{1});

%         tic
        bigvid = zeros([vidObj{1}.Height vidObj{1}.Width ...
            floor([1000.*(length(vidObj)-1)+vidObj{end}.Duration.*vidObj{end}.FrameRate]./ds)],'uint8');
        currentFrame = 0;
        dscount = 0;
        for i = 1:length(vidObj)
            while currentFrame < 5001 && hasFrame(vidObj{i})
                tmp = readFrame(vidObj{i});
                dscount = dscount+1;
                if mod(dscount,ds)==1
                    currentFrame = currentFrame+1;
                    bigvid(:,:,currentFrame) = tmp;
                end
            end
        end
%         toc 

        Mr2 = double(bigvid(10:end-10,10:end-10,3000:5000)./(0.9));

        dfF = bsxfun(@minus,Mr2,nanmedian(Mr2,3));
%         dfF = bsxfun(@rdivide,dfF,nanstd(dfF,[],3));
        dfF(dfF<-2) = 0;
        dfF = double((dfF-nanmin(dfF(:)))./range(dfF(:)));

        sz = size(dfF(:,:,1));
        ns = floor(sz(:,1,1)./5);
        mColor = zeros([size(dfF(:,:,1)) 3]);
        for i = 1:length(doColors(:,1))
            if i < length(doColors(:,1))
                mColor((i-1).*ns+1:(i).*ns,:,:) = repmat(permute(doColors(i,:),[1 3 2]), ...
                    [ns length(dfF(1,:,1))]);
            else
                mColor((i-1).*ns+1:end,:,:) = repmat(permute(doColors(i,:),[1 3 2]), ...
                    [length((i-1).*ns+1:length(dfF(:,1,1))) length(dfF(1,:,1))]);
            end
        end
        
        for i = 1:1:1000
%             tmp = cat(3,256.*dfF(:,:,i)./(0.9),zeros([size(dfF(:,:,1)) 2]));
% % %             [imind,cm] = rgb2ind(uint8(repmat(256.*dfF(:,:,i)./(0.9),[1 1 3])),256); 
            [imind,cm] = rgb2ind(uint8(256.*repmat(dfF(:,:,i)./(0.9),[1 1 3]).*mColor),256); 
%             [imind,cm] = rgb2ind(uint8(tmp),256);
            if i == 1
                checkP([outName,'_dFF.gif']);
                imwrite(imind,cm,[outName,'_dFF.gif'],'gif', 'Loopcount',inf,'DelayTime',0.0333); 
            else
                imwrite(imind,cm,[outName,'_dFF.gif'],'gif','DelayTime',0.0333,'WriteMode','append'); 
            end

            [imind,cm] = rgb2ind(uint8(repmat(Mr2(:,:,i),[1 1 3])),256); 
            if i == 1
                checkP([outName,'_raw.gif']);
                imwrite(imind,cm,[outName,'_raw.gif'],'gif', 'Loopcount',inf,'DelayTime',0.0333); 
            else
                imwrite(imind,cm,[outName,'_raw.gif'],'gif','DelayTime',0.0333,'WriteMode','append'); 
            end
        end
        fprintf(['\tDone.']);
    end
end