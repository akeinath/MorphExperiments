function contVecAnalysis(paths)
    
    clc
    close all
    drawnow
    
    
    pause_thresh = 2;
    pvWindow = 60; % Frames for chunking vec for doVectorSim

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
    
    envSize = [17 17];
    pThresh = 0.025;
    envLabel = [{'sq1'} {'sq2'} {'sq3'} ...
        {'g3'} {'g2'} {'g1'}];
    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        doAl = help_getAlignmentID(s.alignment,length(isM),paths(isM));
        alignMap = s.alignment(doAl).alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        amfr = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        COMs = repmat({[]},[1 length(sessions)]);
        SHCs = repmat({[]},[1 length(sessions)]);
        SICs = repmat({[]},[1 length(sessions)]);
        PFSs = repmat({[]},[1 length(sessions)]);
        MFRs = repmat({[]},[1 length(sessions)]);
        simVecs = repmat({[]},[1 length(sessions)]);
        aos = [];
        envs = [];
        tic
        
        slashInds = find(ismember(paths{1},'/'));
        root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
        slashInds = find(ismember(upiece{mi},'/'));
        root = [root '/' upiece{mi}(slashInds(end)+1:end) '/ContinuousAnalyses'];
        
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','exclude');
            slashInds = find(ismember(sessions{si},'/'));
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            amfr{si} = nanmean(gT,2);
            v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            
%             PFSs{si} = permute(nansum(nansum(tm>0.25,1),2),[3 1 2]) ./ ...
%                 nansum(nansum(~isnan(m(:,:,1))));
            isPC{si} = s.processed.splithalf.wholemap_unmatched.p <= pThresh;
            if isfield(s.processed,'exclude')
                isPC{si} = isPC{si}(s.processed.exclude.SFPs,:);
            end
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
            simVecs{si} = chunkVec(gT,pvWindow);
        end  
        durat = toc;
        fprintf([num2str(durat) ' s\n']);
        
        minVecTimes = nanmin(cellfun(@size,simVecs,repmat({2},size(simVecs))));
        uvecSim = nan([length(alignMap{1}(:,1)) minVecTimes length(sessions)]);
        for si = 1:length(sessions)
            uvecSim(alignMap{1}(:,si)~=0,:,si) = ...
                simVecs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),1:minVecTimes);
        end  
        
        biVecSim = double(uvecSim);
%             biVecSim(isnan(uvecSim)) = 0;

        aIsGood = nan(length(sessions));
        for si = 1:length(sessions) 
            for sj = si+1:length(sessions) 
                a = biVecSim(:,:,si);
                b = biVecSim(:,:,sj);
                aIsGood(si,sj) = nansum(~any(isnan(a),2) & ~any(isnan(b),2));
            end
        end
        doCells = nanmin(aIsGood(:));

        features = [];
        for i = 1:length(envLabel)
            envI = find(ismember(envs,envLabel(i)));
            for k = 1:length(envI)
                features(envI(k),:) = [k i];
            end
        end
        attractorMat = abs(bsxfun(@minus,features(:,2)>=4,[features(:,2)>=4]'));

        mouse = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end);

        if ismember(mouse,{'CAMCA130'})
            morphPoint = [4 4 5 5 6 6];
            vec = [];
            for gi = 1:6
                vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
            end
            attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
                [features(:,2)'>=vec(1:length(features(:,2)))]'));
        end
        
        
        fprintf(['\t\tVector Similarity (independent of place)... '])
        nsims = 30;
        abmMat = nan(length(sessions),length(sessions),nsims);
        abmv = [];
        atxp = [];
        atxc = [];
        admat = [];
        abcMat = nan(length(sessions),length(sessions),nsims);
        apxenv = [];
        strLength = 0;
        groupedVecSim = [];
        for simIter = 1:nsims
            fprintf(repmat('\b',[1 strLength]));
            str = sprintf([ num2str(simIter) ' of ' num2str(nsims)]);
            fprintf(str);
            strLength = length(str);

% % %                 biVecSim = biVecSim(:,:,randperm(length(biVecSim(1,1,:))));

            dmat = nan([length(biVecSim(1,:,1)) length(sessions) length(sessions)]);
            for si = 1:length(sessions) 
                for sj = si+1:length(sessions) 
                    a = biVecSim(:,:,si);
                    b = biVecSim(:,:,sj);
                    isGood = ~any(isnan(a),2) & ~any(isnan(b),2);

                    ac = find(isGood);
                    ac = ac(randperm(length(ac)));
                    isGood = false(size(isGood));
                    isGood(ac(1:doCells)) = true; % Downsample to match number of cells

                    xc = corr(a(isGood,:),b(isGood,:));

%                         xc = squareform(pdist([a(isGood,:) b(isGood,:)]','cosine'));
%                         xc = xc(1:length(a(1,:)),length(a(1,:))+1:end)';

                    dmat(:,sj,si) = nanmax(xc,[],2);
                    dmat(:,si,sj) = nanmax(xc,[],1); % nanmin(nanmean(abs(bsxfun(@minus,a(isGood),b(isGood,:))),1));
                end
            end
            
            
%             tmp = [{[]} {[]}];
%             for i = 1:length(attractorMat)
%                 tmp{1} = [tmp{1}; nanmax(dmat(:,~logical(attractorMat(:,i)),i),[],2)];
%                 tmp{2} = [tmp{2}; nanmax(dmat(:,logical(attractorMat(:,i)),i),[],2)];
%             end
%             groupedVecSim = [groupedVecSim; tmp];
%             mkWhisker(tmp,[{'In'} {'Out'}]);
            
            [a bestSession] = nanmax(dmat,[],2);
            bms = permute(bestSession,[1 3 2]);
            admat = cat(4,admat,dmat);


            meanMaxCorr = nan(length(sessions));
            for i = 1:length(sessions)
                for j = 1:length(sessions)
                    meanMaxCorr(i,j) = nanmean(dmat(:,i,j));
                end
            end
            abcMat(:,:,simIter) = meanMaxCorr;

            timeByCorr = nan(length(sessions).*2-1,length(sessions));
            for i = 1:length(sessions)
                timeByCorr(length(sessions)-i+1:2.*length(sessions)-i,i) = ...
                    meanMaxCorr(:,i);
            end
            atxc = [atxc; nanmean(timeByCorr')];

            bmv = permute(a,[1 3 2]);

            bmvMat = nan(1,length(sessions));
            bmMat = nan(length(sessions));
            for i = 1:length(sessions)
                bmMat(i,:) = (nanmean(bms==i,1));
                bmvMat(i) = (nanmean(bmv(:,i),1));
            end
            bmMat(logical(eye(size(bmMat)))) = nan;
            abmMat(:,:,simIter) = bmMat;
            abmv = [abmv; bmvMat];

            timeByPerc = nan(length(sessions).*2-1,length(sessions));
            for i = 1:length(sessions)
                timeByPerc(length(sessions)-i+1:2.*length(sessions)-i,i) = ...
                    abmMat(:,i,simIter);
            end
            atxp = [atxp; nanmean(timeByPerc')];


            [a condGroup] = ismember(envs,envLabel);
            tmp = repmat({[]},[nanmax(condGroup) nanmax(condGroup)]);
            for i = 1:length(sessions)
                for j = 1:nanmax(condGroup)
                    tmp{condGroup(i),j} = [tmp{condGroup(i),j}; ...
                        bmMat(condGroup==j,i)];
                end
            end
            tmp = cellfun(@nanmean,tmp);
            apxenv = cat(3,apxenv,tmp);
        end

%         tmp = [];
%         for i = 1:2
%             tmp{i} = nanmedian(cat(2,groupedVecSim{:,i}),2);
%             tmp{i} = tmp{i}(length(tmp{i})./2+1:end);
%         end
%         mkWhisker(tmp,[{'In'} {'Out'}]);
        
        close all
        figure
        set(gcf,'position',[50 50 900 250])
        subplot(1,3,1)
        mkGraph(atxp,[],[0.3 0.5 1]);
        hold on
        set(gca,'xtick',[2:5:length(sessions).*2-1], ...
            'xticklabel',[2:5:length(sessions).*2-1]-length(sessions))
        plot([length(sessions) length(sessions)],get(gca,'ylim'),'color','k','linestyle','--')
        xlabel('Lag (days)')
        ylabel('Likelihood of Best Vector Match')
        subplot(1,3,2)
        imagesc(nanmean(apxenv,3));
        hold off
        drawnow
        subplot(1,3,3)
        plot(nanmean(abmv,1));
        hold off
        drawnow


        saveFig(gcf,[root '/VectorSimilarityAnalysis' ...
            num2str(pvWindow./30) 's'],[{'pdf'} {'tiff'}]); % num2str(pvWindow./30) 

        tmp = nanmean(abmMat,3);

        toPlot = [];
        tam = ~logical(attractorMat);
        for i = 1:length(tmp(1,:))
            toPlot = [toPlot; nansum(tmp(tam(:,i),i)) nansum(tmp(~tam(:,i),i))];
        end
        figure
        set(gcf,'position',[50 50 200 300])
        mkWhisker(toPlot,[{'In'} {'Out'}]);
        hold on
        set(gca,'ylim',get(gca,'ylim')+[-0.025 0.025])
        plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--');
        [pval h stats] = signrank(toPlot(:,1),0.5);
        xlabel('Env. Attractor Group')
        ylabel('Likelihood of Best Vector Match')
        set(gca,'ylim',[0 1])
        textY = get(gca,'ylim');
        textY = textY(1)+(diff(textY).*0.05);


        tmp = nanmean(abcMat,3);
        toPlot = [];
        tam = ~logical(attractorMat);
        for i = 1:length(tmp(1,:))
            toPlot = [toPlot; nanmean(tmp(tam(:,i),i)) nanmean(tmp(~tam(:,i),i))];
        end
        figure
        set(gcf,'position',[50 50 200 300])
        mkWhisker(toPlot,[{'In'} {'Out'}]);
        hold on
        set(gca,'ylim',get(gca,'ylim')+[-0.025 0.025])
%             plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--');
        [pval blah stats] = signrank(toPlot(:,1),toPlot(:,2));
        xlabel('Env. Attractor Group')
        ylabel('Likelihood of Best Vector Match')
        textY = get(gca,'ylim');
        textY = textY(1)+(diff(textY).*0.05);
%             tStr = sprintf('Z = %0.2f, p = %0.3f',[stats.zval pval]);
%             text(0.6,textY,tStr)

        saveFig(gcf,[root '/VectorSimilarityAnalysis_MaxCorr_ByAttractor' ...
            num2str(pvWindow./30) 's'],[{'pdf'} {'tiff'}]); % num2str(pvWindow./30) 
end




















