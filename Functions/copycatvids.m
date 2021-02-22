function copycatvids(inFolder,outFolder,prefix)
    if nargin < 3 || isempty(prefix)
        prefix = 'behavCam';
    end
    
    
    paths = getFilePaths(inFolder,'.avi');
    unFolder = outFolder;

    clc
    
    fprintf('Concatenating and copying unpacked files:\n');
    
    tmp = [];
    for i = length(paths):-1:1
        pind = find(ismember(paths{i},'/'),1,'last');
        if length(paths{i}(pind+1:end))>=length(prefix) && ...
                ismember({paths{i}(pind+1:pind+length(prefix))},prefix)
            tmp = [tmp; {paths{i}(1:pind-1)}];
        else
            paths(i) = [];
        end
    end
    tmp = flipud(tmp);

    [a b] = unique(tmp);
    [a groupID] = ismember(tmp,a');
    
    for i = 1:nanmax(groupID)
        doS = paths(i==groupID);
        vidN = nan(length(doS),1);
        for j = 1:length(doS)
            ind = strfind(doS{j},prefix);
            vidN(j) = str2num(doS{j}(ind+8:find(ismember(doS{j},'.'),1,'last')-1));
        end
        [a b] = sort(vidN);
        doS = doS(b);
        
        op = doS{j};
        slashInds = find(ismember(op,'/'));
        tmp = op(1:slashInds(end)-1);
%         tmp(ismember(tmp,'/')) = '_';
        np = [unFolder '/' tmp '.avi'];
        fprintf(['\n\t' tmp]);
        checkP(np);
        
        tic
        writerObj = VideoWriter(np);
        open(writerObj)
        for j = 1:length(doS)
            v = VideoReader(doS{j});    
            frame = 1;
            nF = v.NumberOfFrames;
            v = VideoReader(doS{j});  
            vid = VideoReader(doS{j});    
            while hasFrame(v) && frame<=nF
                thisFrame = read(vid, frame);
                writeVideo(writerObj,thisFrame);
                frame = frame+1;
            end
        end
        close(writerObj)
        tocString = toc;
        fprintf(['\n\t\tTime to complete:  ' num2str(tocString) '\n'])
    end
    
    fprintf('\n\nComplete.\n\n')
    
end