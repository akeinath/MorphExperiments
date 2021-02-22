function alignNwiseSessions(paths,nFold)
    if nargin < 2 || isempty(nFold)
        nFold = 2;
    end
    
    iters = 5;
    
    warning ('off','all');
    
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
    
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);      
        
        combs = nchoosek(1:length(sessions),nFold);

%         %%% Sliding temporal window
%         combs = repmat(1:length(sessions),[length(sessions) 1]);
%         for i = 1:length(sessions)
%             combs(i,:) = circshift(combs(i,:),[0 -(i-1)]);
%         end
%         combs(end-nFold+2:end,:) = [];
%         combs = combs(:,1:nFold);
       
        allPrepped = repmat({[]},[1 length(sessions)]);
        for si = 1:length(sessions)
            fprintf(['\t\tSession:  ' sessions{si}(find(ismember(sessions{si},'/'),1,'last')+1:end) '\n'])
            ref = load(sessions{si},'calcium','processed');       
            prepped = ref.calcium.SFPs .* ...
                    bsxfun(@gt,ref.calcium.SFPs,0.5.*nanmax(nanmax(ref.calcium.SFPs,[],1),[],2));
%             prepped = ref.calcium.SFPs.^(2);
            allPrepped{si} = permute(prepped,[3 1 2]); %msExtractSFPs(ref.calcium);
            if isfield(ref.processed,'exclude')
                allPrepped{si} = allPrepped{si}(ref.processed.exclude.SFPs,:,:);
            end
            if si == 1
                ref = load(sessions{si},'alignment');
                if nFold == 2
                    oldScores = repmat({[]},[length(isM) length(isM)]);
                    oldAlignmentMap = repmat({[]},[length(isM) length(isM)]);
                else
                    oldScores = repmat({[]},[length(combs(:,1)) 1]);
                    oldAlignmentMap = repmat({[]},[length(combs(:,1)) 1]);
                end
                if isfield(ref,'alignment')
                    alignID = help_getAlignmentID(ref.alignment,nFold,paths(isM));
                    if ~isnan(alignID)
                        oldAlignmentMap = ref.alignment(alignID).alignmentMap;
                        if isfield(ref.alignment(alignID),'scores')
                            oldScores = ref.alignment(alignID).scores;
                        end
                    end
                end
            end
        end
        
        tas = repmat({[]},[length(combs(:,1)) 1]);
        mspiece = spiece(isM);
        alignmentMap = repmat({[]},[length(combs(:,1)) 1]);
        for i = 1:length(combs(:,1))
            itermap = repmat({[]},[1 iters]);
            iterstruct = repmat({[]},[1 iters]);
            for si = combs(i,:)              
                prepped = allPrepped{si};
                outP = ['Segments/' upiece{mi} '/' mspiece{si}];% num2str(find(si == combs(i,:)))];
                checkP(outP)
                save(outP,'prepped'); 
            end

            for iteration = 1:iters
                try
                    [map regStruct] = registerCells(['Segments/' upiece{mi} ]);
                catch
                    map = [];
                    regStruct = [];
                end
                itermap{iteration} = map;
                iterstruct{iteration} = regStruct;
                close all
                close all hidden
                drawnow
            end
            try
                rmdir(['Segments/' upiece{mi} ],'s');
            end
            iterstruct = iterstruct(~cellfun(@isempty,itermap));
            itermap = itermap(~cellfun(@isempty,itermap));
            
            scores = [];
            for j = 1:length(iterstruct)
%                 if nFold == 2
%                     scores(j) = nansum(iterstruct{j}.cell_scores( ...
%                         all(iterstruct{j}.cell_to_index_map~=0,2))>0.95);
%                 end
                scores(j) = nanmean(iterstruct{j}.cell_scores);
            end
            
            [a b] = nanmin(cellfun(@length,itermap));
%             [a b] = nanmax(scores);
            ts = [];
            if ~isempty(b)
                map = itermap{b};
                ts = iterstruct{b}.cell_scores;
            end
            alignmentMap{i} = map;
            tas{i} = ts;
        end
        
        sm = repmat({[]},length(sessions));
        nm = repmat({[]},length(sessions));
        for i = 1:length(combs(:,1))
            if nFold == 2
                doAddMap = [alignmentMap(i) oldAlignmentMap(combs(i,1),combs(i,2))];
                doAddScores = [tas(i) oldScores(combs(i,1),combs(i,2))];

                best = [{tas{i}} ...
                    {oldScores{combs(i,1),combs(i,2)}}];
            else
                doAddMap = [alignmentMap(i) oldAlignmentMap(i)];
                doAddScores = [tas(i) oldScores(i)];

                best = [{tas{i}} {oldScores{i}}];
            end
%             best = cellfun(@nanmean,cellfun(@gt,best,[{[0.95]} {[0.95]}],'uniformoutput',false));
%             best = cellfun(@nanmean,best);
            best = cellfun(@length,best);
            
            
%             if nFold == 2
%                 doAddMap = [alignmentMap(i) oldAlignmentMap(combs(i,1),combs(i,2))];
%                 doAddScores = [tas(i) oldScores(combs(i,1),combs(i,2))];
% 
%                 best = [{tas{i}(all(doAddMap{1}~=0,2),:)} ...
%                     {oldScores{combs(i,1),combs(i,2)}(all(doAddMap{2}~=0,2),:)}];
%             else
%                 doAddMap = [alignmentMap(i) oldAlignmentMap(i)];
%                 doAddScores = [tas(i) oldScores(i)];
% 
%                 best = [{tas{i}(all(doAddMap{1}~=0,2),:)} ...
%                     {oldScores{i}(all(doAddMap{2}~=0,2),:)}];
%             end
%             best = cellfun(@nansum,cellfun(@gt,best,[{[0.95]} {[0.95]}],'uniformoutput',false));
            if isempty(best)
                continue
            end
            best(best==0) = nan;
            [a b] = nanmax(best);

            if nFold == 2
                nm{combs(i,1),combs(i,2)} = doAddMap{b};
                sm{combs(i,1),combs(i,2)} = doAddScores{b};
            else
                nm{i} = doAddMap{b};
                sm{i} = doAddScores{b};
            end
        end
        alignmentMap = nm;
        scoreMap = sm;

        ref = load(sessions{1});
        if isfield(ref,'alignment')
            doInd = help_getAlignmentID(ref.alignment,nFold,paths);
            if isnan(doInd)
                doInd = length(ref.alignment)+1;
            end
        else
            doInd = 1;
        end
        
        ref.alignment(doInd).alignmentMap = alignmentMap;
        ref.alignment(doInd).scores = scoreMap;
        ref.alignment(doInd).combs = combs;
        ref.alignment(doInd).sessions = sessions;
        ref.alignment(doInd).nFold = nFold;
        save(sessions{1},'-struct','ref','-v7.3');
    end
end