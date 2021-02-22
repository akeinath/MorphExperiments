function [] = quantWaves(folderpath)

    ds = 1;

    clc

    filePrefix = 'msvideo';
    
    ap = getFilePaths(folderpath,'.avi');
    for i = length(ap):-1:1
        slind = find(ismember(ap{i},'\/'));
        if ~ismember({ap{i}(slind(end)+1:end-4)},{filePrefix})
            ap(i) = [];
        end
    end
    
    ap = ap(2);
    
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
            while hasFrame(vidObj{i})
                tmp = readFrame(vidObj{i});
                
                currentFrame = currentFrame+1;
                bigvid(:,:,currentFrame) = tmp;
            end
        end

        bigvid = double(bigvid);

        dfF = bsxfun(@minus,bigvid,nanmedian(bigvid,3));
%         dfF = bsxfun(@rdivide,dfF,nanstd(dfF,[],3));
%         dfF(dfF<-2) = 0;
        dfF = double((dfF-nanmin(dfF(:)))./range(dfF(:)));

     
    end
end