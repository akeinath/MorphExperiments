function meanFrameTranslationTests(paths)

    clc
    close all
    drawnow
    
    warning off all
     %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    labels = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
   
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        
        amf = repmat({[]},[1 length(sessions)]);
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)
            s = load(sessions{si},'processed');
            amf{si} = s.processed.meanFrame;
        end  
        
        tmp = cellfun(@size,amf,'uniformoutput',false);
        tmp = cat(1,tmp{:});
        pad2size = nanmax(tmp);
        
        camf = cat(3,amf{:});
        for i = 1:length(sessions)
            camf(:,:,i) = normMeanFrame(camf(:,:,i),8);
        end
        
        [x y] = meshgrid([1:length(camf(1,:,1))]-length(camf(1,:,1))./2+0.5, ...
        [1:length(camf(:,1,1))]-length(camf(:,1,1))./2+0.5);
        d2c = sqrt(x.^2 + y.^2);
        mask = d2c < nanmin(size(camf(:,:,1))).*0.4;

        for i = 1:length(sessions)
            for j = i+1:length(sessions)
                ta = camf(:,:,i);
                tb = camf(:,:,j);
                ta(~mask) = nan;
                tb(~mask) = nan;
                xc = pvxcorr(ta,tb,[25 25]);
            end            
        end

        title('filtered gray scale image');
    end
 end    