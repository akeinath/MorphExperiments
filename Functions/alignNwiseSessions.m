function alignNwiseSessions(paths,nFold)
    if nargin < 2 | isempty(nFold)
        nFold = 2;
    end
    
    warning ('off','all');
    
    
    iters = 5;
    
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
        
        if nFold > length(sessions)
            continue
        end
        
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
                    oldAlignmentMap = repmat({[]},[length(isM) length(isM)]);
                else
                    oldAlignmentMap = repmat({[]},[length(combs(:,1))]);
                end
                if isfield(ref,'alignment')
                    alignID = help_getAlignmentID(ref.alignment,nFold,paths(isM));
                    if ~isnan(alignID)
                        oldAlignmentMap = ref.alignment(alignID).alignmentMap; % Comment out to NOT combine
                    end
                end
            end
        end
        
        mspiece = spiece(isM);
        alignmentMap = repmat({[]},[length(combs(:,1)) 1]);
        for i = 1:length(combs(:,1))
            
%             if nFold==2 
%                 if ~isempty(oldAlignmentMap{combs(i,1),combs(i,2)})
%                     continue
%                 end
%             else
%                 if ~isempty(oldAlignmentMap{i})
%                     continue
%                 end
%             end
            
            itermap = repmat({[]},[1 iters]);
            for si = combs(i,:)              
                prepped = allPrepped{si};
                outP = ['SegmentsForAlignment/' mspiece{si}];% num2str(find(si == combs(i,:)))];
                checkP(outP)
                save(outP,'prepped'); 
            end

            for iteration = 1:iters
                try
                    [map regStruct] = registerCells(['SegmentsForAlignment']);
                catch
                    map = [];
                end
                itermap{iteration} = map;
                close all
                close all hidden
                drawnow
            end
            rmdir(['SegmentsForAlignment'],'s');
            
            itermap = itermap(~cellfun(@isempty,itermap));
            
            if ~isempty(itermap)
                [a b] = nanmin(cellfun(@length,itermap));
                map = itermap{b};
                alignmentMap{i} = map;
            end
        end
        if nFold == 2
            nm = repmat({[]},length(sessions));
            for i = 1:length(combs(:,1))
                best = [alignmentMap(i) oldAlignmentMap(combs(i,1),combs(i,2))];
                best = best(~cellfun(@isempty,best));
                if isempty(best)
                    continue
                end
                [a b] = nanmin(cellfun(@length,best));
                nm{combs(i,1),combs(i,2)} = best{b};
            end
            alignmentMap = nm;
        else
            nm = repmat({[]},[length(combs(:,1)) 1]);
            for i = 1:length(combs(:,1))
                best = [alignmentMap(i) oldAlignmentMap(i)];
                best = best(~cellfun(@isempty,best));
                if isempty(best)
                    continue
                end
                [a b] = nanmin(cellfun(@length,best));
                nm{i} = best{b};
            end
            alignmentMap = nm;
        end
        
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
        ref.alignment(doInd).combs = combs;
        ref.alignment(doInd).sessions = sessions;
        ref.alignment(doInd).nFold = nFold;
        if nFold == length(sessions)
            ref.alignment(doInd).regDetails = regStruct;
        end
        save(sessions{1},'-struct','ref','-v7.3');
    end
end