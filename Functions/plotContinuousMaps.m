function plotContinuousMaps(paths)

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
    
    envSize = [17 19];
    
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        doAl = help_getAlignmentID(s.alignment,length(isM),paths(isM));
        alignMap = s.alignment(doAl).alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        envs = [];
        tic
        
        slashInds = find(ismember(paths{1},'/'));
        root = ['Plots/' paths{1}(slashInds(1):slashInds(2)-1)];
        slashInds = find(ismember(upiece{mi},'/'));
        root = [root '/' upiece{mi}(slashInds(end)+1:end) '/ContinuousAnalyses'];
        close all
        drawnow;
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)

            s = load(sessions{si},'processed','exclude');

            slashInds = find(ismember(sessions{si},'/'));
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            [m os unm] = mkTraceMaps(s.processed.p,gT,[],envSize);
            am{si} = m;
            
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
        end  

        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        
        for si = 1:length(sessions)
            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
        end
        
        catMaps = nan(length(um(:,1,1,1)),length(um(1,1,1,:)).*length(um(1,:,1,1)),length(um(1,1,:,1)));
        for si = 1:length(sessions)
            for k = 1:length(um(1,1,:,1))
                catMaps(:,(si-1).*length(um(1,:,1,1))+1:(si).*length(um(1,:,1,1)),k) = ...
                    um(:,:,k,si)./nanmax(nanmax(um(:,:,k,si)));
            end
        end
        
        doSize = [6 2];
        for gi = 1:floor(length(catMaps(1,1,:))./prod(doSize))
            figure
            set(gcf,'position',[50 50 fliplr(doSize).*[200 40].*4])
            for k = 1:prod(doSize)
                try
                    dm = catMaps(:,:,(gi-1).*prod(doSize)+k);
                    subplot(doSize(1),doSize(2),k)
                    imagesc(dm)
                    axis equal
                    axis off
                    alpha(double(~isnan(dm)))
                end
            end
            saveFig(gcf,[root '/ContinuousMaps/Part_' num2str(gi)],[{'pdf'} {'tiff'}])
        end
    end
end