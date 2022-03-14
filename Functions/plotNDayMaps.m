function plotNDayMaps(paths)
    
    warning off all
    
    clc
    fprintf('\n')
    %%% Reliability constaint
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    envSize = [17 17];
    doPlot = true;
    pause_thresh = -2;
    envLabel = [{'sq1'} {'sq2'} {'sq3'} ...
        {'g3'} {'g2'} {'g1'}];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)},'alignment');
        doAl = help_getAlignmentID(s.alignment,length(isM),paths(isM));
        alignMap = s.alignment(doAl).alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        
        envs = [];
        tic
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','exclude');
            slashInds = find(ismember(sessions{si},'/'));
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            [m os] = mkTraceMaps(s.processed.p,gT,v>=pause_thresh,envSize);
            am{si} = m;
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
        end  
        durat = toc;
        fprintf([num2str(durat) ' s\n']);
        
        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        for si = 1:length(sessions)
            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
        end  
        um = um(:,:,all(~(all(all(isnan(um),1),2)),4),:);
        
        [a b] = ismember(envs,envLabel);
        [a order] = sort(b);
        
        outP = ['Plots/NWiseAlignedCellMaps/' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '/' ...
            'AlignedMaps'];
        helpPlotMaps(um(:,:,:,[order]),outP);

    end
end