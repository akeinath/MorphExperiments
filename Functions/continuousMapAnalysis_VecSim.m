function continuousMapAnalysis(paths,doPlot)
    
    clc
    close all
    drawnow
    
    doMDS = true;
    
    pause_thresh = 2;
    cell_thresh_for_burst = 0.1;
    doPreplay = false;
    doSubBins = false;
    subBinSize = 150;
    doVectorSim = false;
    doAnatomy = false;
    pvWindow = 60; % Frames for chunking vec for doVectorSim
    doCOM = false;
    if nargin < 2
        doPlot = false;
    end

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
    
    allR2Penalty = [];
    allTR2 = [];
    allCellProps = [];

    envSize = [17 17];
    includeLags = inf;
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
        SFPs = s.alignment(doAl).regDetails.spatial_footprints_corrected;
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
            if doAnatomy
                s = load(sessions{si},'processed','exclude','calcium');
            else
                s = load(sessions{si},'processed','exclude');
            end
            slashInds = find(ismember(sessions{si},'/'));
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            amfr{si} = nanmean(gT,2);
            v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            [m os] = mkTraceMaps(s.processed.p,gT,v>=pause_thresh,envSize);
            aos = cat(3,aos,os);
            am{si} = m;
            
            tm = m./repmat(nanmax(nanmax(m,[],1),[],2),size(m(:,:,1)));
%             PFSs{si} = permute(nansum(nansum(tm>0.25,1),2),[3 1 2]) ./ ...
%                 nansum(nansum(~isnan(m(:,:,1))));
            PFSs{si} = permute(nansum(nansum(tm>0.25,1),2),[3 1 2]);
            SICs{si} = s.processed.splithalf.wholemap_unmatched_si.val;
            SHCs{si} = s.processed.splithalf.wholemap_unmatched.val;
            isPC{si} = s.processed.splithalf.wholemap_unmatched.p <= pThresh;
            MFRs{si} = nanmean(gT(:,v>=pause_thresh),2);
            if isfield(s.processed,'exclude')
                isPC{si} = isPC{si}(s.processed.exclude.SFPs,:);
            end
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
            if doVectorSim    
                simVecs{si} = chunkVec(gT,pvWindow);
            end
            if doPreplay
                [preVecs{si} pvals] = getBurstVectors(s.processed.p,gT, ...
                    pause_thresh,cell_thresh_for_burst,isPC{si});
            end
            if doCOM
                COMs{si} = getCOMs(m);
            end
            if doAnatomy
                MeanFrames{si} = s.calcium.meanFrame;
                SFPs{si} = s.calcium.SFPs(:,:,s.processed.exclude.SFPs);
            end
        end  
        durat = toc;
        fprintf([num2str(durat) ' s\n']);
        
        minVecTimes = nanmin(cellfun(@size,simVecs,repmat({2},size(simVecs))));
        uvecSim = nan([length(alignMap{1}(:,1)) minVecTimes length(sessions)]);
        umfr = nan([length(alignMap{1}(:,1)) length(sessions)]);
        upfs = nan([length(alignMap{1}(:,1)) length(sessions)]);
        usic = nan([length(alignMap{1}(:,1)) length(sessions)]);
        ushc = nan([length(alignMap{1}(:,1)) length(sessions)]);
        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        
        tmp = size(SFPs{1}(1,:,:));
        if doAnatomy
            usfps = nan([tmp(2:3) length(alignMap{1}(:,1)) length(sessions)]);
        end
        for si = 1:length(sessions)
%             isCPC = isPC{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
%             tmp = false(length(alignMap{1}(:,si)),1);
%             tmp(alignMap{1}(:,si)~=0) = isCPC;
%             alignMap{1}(~tmp,si) = 0; %%% Remove the non-place cells
            
            if doVectorSim
                uvecSim(alignMap{1}(:,si)~=0,:,si) = ...
                    simVecs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),1:minVecTimes);
            end
            
            umfr(alignMap{1}(:,si)~=0,si) = MFRs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            upfs(alignMap{1}(:,si)~=0,si) = PFSs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            usic(alignMap{1}(:,si)~=0,si) = SICs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            ushc(alignMap{1}(:,si)~=0,si) = SHCs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            
            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
            
            if doAnatomy
                SFPs{si} = permute(SFPs{si},[2 3 1]);
                usfps(:,:,alignMap{1}(:,si)~=0,si) = SFPs{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
            end
            if doPreplay
                tmp = nan(length(alignMap{1}(:,1)),length(preVecs{si}(1,:)));
                tmp(alignMap{1}(:,si)~=0,:) = preVecs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),:);
                preVecs{si} = tmp;
            end
            if doCOM
                tmp = nan(length(alignMap{1}(:,1)),length(COMs{si}(1,:)));
                tmp(alignMap{1}(:,si)~=0,:) = COMs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),:);
                COMs{si} = tmp;
            end
        end
        
        if doAnatomy
            msfps = nanmean(usfps,4);
        end
        
% % %         if doVectorSim
% % %             biVecSim = double(uvecSim);
% % % %             biVecSim(isnan(uvecSim)) = 0;
% % %                
% % %             aIsGood = nan(length(sessions));
% % %             for si = 1:length(sessions) 
% % %                 for sj = si+1:length(sessions) 
% % %                     a = biVecSim(:,:,si);
% % %                     b = biVecSim(:,:,sj);
% % %                     aIsGood(si,sj) = nansum(~any(isnan(a),2) & ~any(isnan(b),2));
% % %                 end
% % %             end
% % %             doCells = nanmin(aIsGood(:));
% % % %             doCells = 31;
% % % 
% % %             fprintf(['\t\tVector Similarity (independent of place)... '])
% % %             nsims = 30;
% % %             abmMat = nan(length(sessions),length(sessions),nsims);
% % %             abmv = [];
% % %             atxp = [];
% % %             atxc = [];
% % %             admat = [];
% % %             abcMat = nan(length(sessions),length(sessions),nsims);
% % %             apxenv = [];
% % %             strLength = 0;
% % %             for simIter = 1:nsims
% % %                 fprintf(repmat('\b',[1 strLength]));
% % %                 str = sprintf([ num2str(simIter) ' of ' num2str(nsims)]);
% % %                 fprintf(str);
% % %                 strLength = length(str);
% % %                 
% % % % % %                 biVecSim = biVecSim(:,:,randperm(length(biVecSim(1,1,:))));
% % %                 
% % %                 dmat = nan([length(biVecSim(1,:,1)) length(sessions) length(sessions)]);
% % %                 for si = 1:length(sessions) 
% % %                     for sj = si+1:length(sessions) 
% % %                         a = biVecSim(:,:,si);
% % %                         b = biVecSim(:,:,sj);
% % %                         isGood = ~any(isnan(a),2) & ~any(isnan(b),2);
% % % 
% % %                         ac = find(isGood);
% % %                         ac = ac(randperm(length(ac)));
% % %                         isGood = false(size(isGood));
% % %                         isGood(ac(1:doCells)) = true; % Downsample to match number of cells
% % % 
% % %                         xc = corr(a(isGood,:),b(isGood,:));
% % %                         
% % % %                         xc = squareform(pdist([a(isGood,:) b(isGood,:)]','cosine'));
% % % %                         xc = xc(1:length(a(1,:)),length(a(1,:))+1:end)';
% % %                         
% % %                         dmat(:,sj,si) = nanmax(xc,[],2);
% % %                         dmat(:,si,sj) = nanmax(xc,[],1); % nanmin(nanmean(abs(bsxfun(@minus,a(isGood),b(isGood,:))),1));
% % %                     end
% % %                 end
% % %                 [a bestSession] = nanmax(dmat,[],2);
% % %                 bms = permute(bestSession,[1 3 2]);
% % %                 admat = cat(4,admat,dmat);
% % %                 
% % %                 
% % %                 meanMaxCorr = nan(length(sessions));
% % %                 for i = 1:length(sessions)
% % %                     for j = 1:length(sessions)
% % %                         meanMaxCorr(i,j) = nanmean(dmat(:,i,j));
% % %                     end
% % %                 end
% % %                 abcMat(:,:,simIter) = meanMaxCorr;
% % %                 
% % %                 timeByCorr = nan(length(sessions).*2-1,length(sessions));
% % %                 for i = 1:length(sessions)
% % %                     timeByCorr(length(sessions)-i+1:2.*length(sessions)-i,i) = ...
% % %                         meanMaxCorr(:,i);
% % %                 end
% % %                 atxc = [atxc; nanmean(timeByCorr')];
% % %                 
% % %                 bmv = permute(a,[1 3 2]);
% % %                 
% % %                 bmvMat = nan(1,length(sessions));
% % %                 bmMat = nan(length(sessions));
% % %                 for i = 1:length(sessions)
% % %                     bmMat(i,:) = (nanmean(bms==i,1));
% % %                     bmvMat(i) = (nanmean(bmv(:,i),1));
% % %                 end
% % %                 bmMat(logical(eye(size(bmMat)))) = nan;
% % %                 abmMat(:,:,simIter) = bmMat;
% % %                 abmv = [abmv; bmvMat];
% % %                 
% % %                 timeByPerc = nan(length(sessions).*2-1,length(sessions));
% % %                 for i = 1:length(sessions)
% % %                     timeByPerc(length(sessions)-i+1:2.*length(sessions)-i,i) = ...
% % %                         abmMat(:,i,simIter);
% % %                 end
% % %                 atxp = [atxp; nanmean(timeByPerc')];
% % %                              
% % %                 
% % % % % %                 [a condGroup] = ismember(envs,envLabel);
% % % % % %                 minSamps = nanmin(histc(condGroup,[1:nanmax(condGroup)]))-1;
% % % % % %                 tmp = repmat({[]},[nanmax(condGroup) nanmax(condGroup)]);
% % % % % %                 for i = 1:length(sessions)
% % % % % %                     
% % % % % %                     doInclude = false(length(condGroup),1);
% % % % % %                     for j = 1:nanmax(condGroup)
% % % % % %                         ginds = find(condGroup==j);
% % % % % %                         ginds = ginds(randperm(length(ginds)));
% % % % % %                         while j==condGroup(i) && ismember(i,ginds(1:minSamps))
% % % % % %                             ginds = ginds(randperm(length(ginds)));
% % % % % %                         end
% % % % % %                         doInclude(ginds(1:minSamps)) = true;
% % % % % %                     end
% % % % % %                     
% % % % % %                     tdmat = dmat;
% % % % % %                     tdmat(:,~doInclude,:) = nan;
% % % % % %                     
% % % % % %                     [a bestSession] = nanmax(tdmat,[],2);
% % % % % %                     tbms = permute(bestSession,[1 3 2]);
% % % % % %                     tbmMat = nan(length(sessions));
% % % % % %                     for j = 1:length(sessions)
% % % % % %                         tbmMat(j,:) = (nanmean(tbms==j,1));
% % % % % %                     end
% % % % % %                     for j = 1:nanmax(condGroup)
% % % % % %                         tmp{condGroup(i),j} = [tmp{condGroup(i),j}; ...
% % % % % %                             nansum(tbmMat(condGroup==j,i))];
% % % % % %                     end
% % % % % %                 end
% % % % % %                 tmp = cellfun(@nanmean,tmp);
% % % % % %                 apxenv = cat(3,apxenv,tmp);
% % %                 
% % %                 
% % %                 [a condGroup] = ismember(envs,envLabel);
% % %                 tmp = repmat({[]},[nanmax(condGroup) nanmax(condGroup)]);
% % %                 for i = 1:length(sessions)
% % %                     for j = 1:nanmax(condGroup)
% % %                         tmp{condGroup(i),j} = [tmp{condGroup(i),j}; ...
% % %                             bmMat(condGroup==j,i)];
% % %                     end
% % %                 end
% % %                 tmp = cellfun(@nanmean,tmp);
% % %                 apxenv = cat(3,apxenv,tmp);
% % %             end
% % %             
% % %             close all
% % %             figure
% % %             set(gcf,'position',[50 50 900 250])
% % %             subplot(1,3,1)
% % %             mkGraph(atxp,[],[0.3 0.5 1]);
% % %             hold on
% % %             set(gca,'xtick',[2:5:length(sessions).*2-1], ...
% % %                 'xticklabel',[2:5:length(sessions).*2-1]-length(sessions))
% % %             plot([length(sessions) length(sessions)],get(gca,'ylim'),'color','k','linestyle','--')
% % %             xlabel('Lag (days)')
% % %             ylabel('Likelihood of Best Vector Match')
% % %             subplot(1,3,2)
% % %             imagesc(nanmean(apxenv,3));
% % %             hold off
% % %             drawnow
% % %             subplot(1,3,3)
% % %             plot(nanmean(abmv,1));
% % %             hold off
% % %             drawnow
% % %             
% % % 
% % %             saveFig(gcf,[root '/VectorSimilarityAnalysis' ...
% % %                 num2str(pvWindow./30) 's'],[{'pdf'} {'tiff'}]); % num2str(pvWindow./30) 
% % %             
% % %             features = [];
% % %             for i = 1:length(envLabel)
% % %                 envI = find(ismember(envs,envLabel(i)));
% % %                 for k = 1:length(envI)
% % %                     features(envI(k),:) = [k i];
% % %                 end
% % %             end
% % %             attractorMat = abs(bsxfun(@minus,features(:,2)>=4,[features(:,2)>=4]'));
% % % 
% % %             mouse = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end);
% % % 
% % %             if ismember(mouse,{'CAMCA130'})
% % %                 morphPoint = [4 4 5 5 6 6];
% % %                 vec = [];
% % %                 for gi = 1:6
% % %                     vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
% % %                 end
% % %                 attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
% % %                     [features(:,2)'>=vec(1:length(features(:,2)))]'));
% % %             end
% % %             
% % %             if ismember(mouse,{'AKCA130'})
% % %                 morphPoint = [4 4 5 5 5 5];
% % %                 vec = [];
% % %                 for gi = 1:6
% % %                     vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
% % %                 end
% % %                 attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
% % %                     [features(:,2)'>=vec(1:length(features(:,2)))]'));
% % %             end
% % %             
% % % %             if ismember(mouse,{'AKCA127'})
% % % %                 morphPoint = [4 4 4 4 4];
% % % %                 vec = [];
% % % %                 for gi = 1:4
% % % %                     vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
% % % %                 end
% % % %                 attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
% % % %                     [features(:,2)'>=vec(1:length(features(:,2)))]'));
% % % %             end
% % %             
% % %             tmp = nanmean(abmMat,3);
% % %             
% % %             toPlot = [];
% % %             tam = ~logical(attractorMat);
% % %             for i = 1:length(tmp(1,:))
% % %                 toPlot = [toPlot; nansum(tmp(tam(:,i),i)) nansum(tmp(~tam(:,i),i))];
% % %             end
% % %             figure
% % %             set(gcf,'position',[50 50 200 300])
% % %             mkWhisker(toPlot,[{'In'} {'Out'}]);
% % %             hold on
% % %             set(gca,'ylim',get(gca,'ylim')+[-0.025 0.025])
% % %             plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--');
% % %             [pval h stats] = signrank(toPlot(:,1),0.5);
% % %             xlabel('Env. Attractor Group')
% % %             ylabel('Likelihood of Best Vector Match')
% % %             set(gca,'ylim',[0 1])
% % %             textY = get(gca,'ylim');
% % %             textY = textY(1)+(diff(textY).*0.05);
% % %             
% % %             
% % % %             tam = ~logical(attractorMat);
% % % %             tamXlag = nan([length(tmp(1,:)) 2]);
% % % %             for lag = 1:length(tmp(1,:))
% % % %                 doLag = diag(diag(true(size(tmp)),lag),lag);
% % % %                 doLag = doLag|doLag';
% % % %                 
% % % %                 tamXlag(lag,:) = [nanmean(tmp(tam&doLag)) nanmean(tmp(~tam&doLag))];
% % % %             end
% % % %             figure
% % % %             plot(tamXlag)
% % % %             drawnow
% % %             
% % % %             admat = permute(admat,[2 3 1 4]);
% % % %             bvmat = nan([length(admat(1,:,1,1)) length(admat(1,1,:,1)) nsims]);
% % % %             for i = 1:nsims
% % % %                 for j = 1:length(admat(1,1,:,1))
% % % %                     for k = 1:length(admat(1,:,1,1))
% % % %                         tmp = admat(:,k,j,i);
% % % %                         bvmat(k,j,i) = nanmax(tmp(tam(:,k)))-nanmax(tmp(~tam(:,k)));
% % % %                     end
% % % %                 end
% % % %             end
% % %             
% % %             
% % %             tmp = nanmean(abcMat,3);
% % %             toPlot = [];
% % %             tam = ~logical(attractorMat);
% % %             for i = 1:length(tmp(1,:))
% % %                 toPlot = [toPlot; nanmean(tmp(tam(:,i),i)) nanmean(tmp(~tam(:,i),i))];
% % %             end
% % %             figure
% % %             set(gcf,'position',[50 50 200 300])
% % %             mkWhisker(toPlot,[{'In'} {'Out'}]);
% % %             hold on
% % %             set(gca,'ylim',get(gca,'ylim')+[-0.025 0.025])
% % % %             plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--');
% % %             [pval blah stats] = signrank(toPlot(:,1),toPlot(:,2));
% % %             xlabel('Env. Attractor Group')
% % %             ylabel('Likelihood of Best Vector Match')
% % %             textY = get(gca,'ylim');
% % %             textY = textY(1)+(diff(textY).*0.05);
% % % %             tStr = sprintf('Z = %0.2f, p = %0.3f',[stats.zval pval]);
% % % %             text(0.6,textY,tStr)
% % %             
% % %             saveFig(gcf,[root '/VectorSimilarityAnalysis_MaxCorr_ByAttractor' ...
% % %                 num2str(pvWindow./30) 's'],[{'pdf'} {'tiff'}]); % num2str(pvWindow./30) 
% % % 
% % % 
% % %             continue
% % %         end
% % %         
% % %         if doPreplay
% % %             tic
% % %             iter = 0;
% % %             strLength = 0;
% % %             fprintf(['\t\tSVM preplay vec classification... '])
% % %             preVecSVMs = nan([length(sessions) length(sessions) length(sessions)]);
% % %             for vecI = 1:length(sessions)
% % %                 for si = 1:length(sessions)  
% % %                     for sj = si+1:length(sessions)   
% % %                         iter = iter+1;
% % %                         fprintf(repmat('\b',[1 strLength]));
% % %                         str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2).*length(sessions))]);
% % %                         fprintf(str);
% % %                         strLength = length(str);
% % %                 
% % %                         
% % %                         cv1 = permute(um(:,:,:,si),[3 1 2]);
% % %                         cv1 = reshape(cv1,[length(cv1(:,1,1)) prod(size(um(:,:,1,1)))]);
% % %                         cv2 = permute(um(:,:,:,sj),[3 1 2]);
% % %                         cv2 = reshape(cv2,[length(cv2(:,1,1)) prod(size(um(:,:,1,1)))]);
% % %                         goodC = ~all(isnan(cv1),2) & ~all(isnan(cv2),2) & ...
% % %                             ~all(isnan(preVecs{vecI}),2);
% % %                         cv1(~goodC,:) = [];
% % %                         cv2(~goodC,:) = [];
% % %                         goodP = ~all(isnan(cv1),1) & ~all(isnan(cv2),1);
% % %                         cv1(:,~goodP) = [];
% % %                         cv2(:,~goodP) = [];
% % %                         
% % %                         tsvm = svmtrain([cv1 cv2]', ...
% % %                             [ones(1,length(cv1(1,:))) 2.*ones(1,length(cv2(1,:)))]');
% % %                         [predLabel] = svmclassify(tsvm,preVecs{vecI}(goodC,:)');
% % %                         preVecSVMs(si,sj,vecI) = nanmean((predLabel-1));
% % %                     end
% % %                 end
% % %             end
% % %             durat = toc;
% % %             fprintf([num2str(durat) ' s\n']);
% % %         end
% % %         
%         if doAnatomy
%             fprintf(['\n\t\tAligning SFPs: '])
%             relShift = nan(length(sessions));
%             for si = round(length(sessions)./2)
%                 a = MeanFrames{si};
%                 fa = fftshift(fft2(a));
%                 [x y] = meshgrid(1:length(a(1,:)),1:length(a(:,1)));
%                 x = x-ceil(length(a(1,:))./2);
%                 y = y-ceil(length(a(:,1))./2);
%                 d2c = sqrt(x.^2 + y.^2);
%                 a2 = ifft2(fftshift(double(d2c>30).*fa));
%                 imagesc(abs(a2))
%                 for sj = si+1:length(sessions)
%                     b = MeanFrames{sj};
%                     fb = fftshift(fft2(b));
%                     [x y] = meshgrid(1:length(b(1,:)),1:length(b(:,1)));
%                     x = x-ceil(length(b(1,:))./2);
%                     y = y-ceil(length(b(:,1))./2);
%                     d2c = sqrt(x.^2 + y.^2);
%                     b2 = ifft2(fftshift(double(d2c>30).*fb));
%                     imagesc(abs(b2))
%                     xc = pvxcorr(abs(a2),abs(b2),[30 30],0);
%                     q = find(xc==nanmax(xc(:)));
%                     [q w] = ind2sub(size(xc),nanmedian(q));
%                     shift = [q w] - round(size(xc)./2);
%                 end
%             end
%         end
        
        kDeformFit = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        kDeformSim = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        sim = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        kCrossCOM = repmat({[]},length(sessions));
        iter = 0;
        strLength = 0;
        fprintf(['\n\t\tComputing pairwise map comparisons: '])
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
                fprintf(str);
                strLength = length(str);
                
                tmp1 = um(:,:,:,si);
                tmp2 = um(:,:,:,sj);
                
% % %                 isOff1 = permute(all(all(isnan(tmp1),1),2),[3 1 2]);
% % %                 isOff2 = permute(all(all(isnan(tmp2),1),2),[3 1 2]);
% % %                 
% % %                 filler1 = nan(size(tmp1(:,:,1)));
% % %                 filler1(~all(isnan(tmp1),3)) = 0;
% % %                 tmp1(:,:,isOff1) = repmat(filler1,[1 1 nansum(isOff1)]);
% % %                 
% % %                 filler2 = nan(size(tmp2(:,:,1)));
% % %                 filler2(~all(isnan(tmp2),3)) = 0;
% % %                 tmp2(:,:,isOff2) = repmat(filler2,[1 1 nansum(isOff2)]);
                
                
                xc = pvxcorr3(tmp1,tmp2,[1 1 1],30);
%                 [a b c] = ind2sub(size(xc),nanmedian(find(xc==nanmax(xc(:)))));
                a = 2;
                b = 2;
                c = 2;
                ivals = xcorr3transform(tmp1,tmp2,[a-ceil(length(xc(:,1,1))./2) ...
                    b-ceil(length(xc(1,:,1))./2) c-ceil(length(xc(1,1,:))./2)]);
                sim(si,sj,:) = ivals;
                
% % %                 isOff1 = permute(all(all(isnan(tmp1),1),2),[3 1 2]);
% % %                 isOff2 = permute(all(all(isnan(tmp2),1),2),[3 1 2]);
% % % 
% % %                 [kDeformSim(si,sj,~isOff1&~isOff2) kDeformFit(si,sj,~isOff1&~isOff2)] = ...
% % %                     fitDeformation(tmp1(:,:,~isOff1&~isOff2), ...
% % %                     tmp2(:,:,~isOff1&~isOff2)); %x y both neither
% % %                 
% % %                 %%% Correct for those that are best fit by no change
% % %                 [a b] = nanmax([permute(kDeformSim(si,sj,~isOff1&~isOff2),[3 1 2]) ...
% % %                     ivals(~isOff1&~isOff2)],[],2);
% % %                 tind = find(~isOff1&~isOff2);
% % %                 kDeformSim(si,sj,tind) = a;
% % %                 kDeformFit(si,sj,tind(b==2)) = 4;
                
                
%                 figure(1)
%                 scatter(permute(sim(si,sj,~isOff1&~isOff2),[3 1 2]), ...
%                     permute(kDeformSim(si,sj,~isOff1&~isOff2),[3 1 2]))
%                 axis equal
%                 axis square
%                 drawnow
                if doCOM
                    kCrossCOM{si,sj} = [COMs{si} COMs{sj}];
                end
                
            end
        end        
        
        percentDeforms = nan([length(sessions) length(sessions) 4]);
        tmp = kDeformFit;
        tmp(kDeformSim<0.7) = nan;
        for i = 1:4
            percentDeforms(:,:,i) = nansum(tmp==i,3)./nansum(~isnan(tmp),3);
        end
        
        percentDeforms = nan([length(sessions) length(sessions) 4]);
        for i = 1:4
            percentDeforms(:,:,i) = nansum(kDeformFit==i,3)./nansum(~isnan(kDeformFit),3);
        end
        
        if doCOM
            [a condGroup] = ismember(envs,envLabel);
            comVecPlot(groupMat(kCrossCOM,condGroup,false),[root '/COM_Vectors']);
        end
        
%         tmp = repmat({[]},[length(envs)]);
%         for i = 1:length(envs)
%             for j = 1:length(envs)
%                 tmp{i,j} = sim(i,j,~isnan(sim(i,j,:)));
%             end
%         end
%         
%         percentDeforms(:,:,4)

        figure
        tmp = squarify(nanmedian(sim,3));
        imagesc(tmp)
        alpha(double(~isnan(tmp)))
%         caxis([0 1])
        colorbar
        axis square
        axis equal
        axis off
        saveFig(gcf,[root '/RDM_All_Pairwise_Sessions'],[{'pdf'} {'tiff'}])
        
%         figure
%         tmp = squarify(nanmedian(kDeformSim,3));
%         imagesc(tmp)
%         alpha(double(~isnan(tmp)))
% %         caxis([0 1])
%         colorbar
%         axis square
%         axis equal
%         axis off
%         saveFig(gcf,[root '/RDM_All_Pairwise_Sessions_GeoCorrected'],[{'pdf'} {'tiff'}])

        if doMDS
            mds2D(sim,envs,envLabel,[root '/MDS_2D']);
%             mds2D(kDeformSim,envs,envLabel,[root '/MDS_GeoCorrected_2D']);
%             mds3D(sim,envs,envLabel,[root '/MDS_3D.gif']);
%             mds3D(kDeformSim,envs,envLabel,[root '/MDS_GeoCorrected_3D.gif']);
            mdsND(sim,envs,envLabel,[root '/MDS_nD/'],5);
%             mdsND(kDeformSim,envs,envLabel,[root '/MDS_nD/GeoCorrected'],5);
        end
        
        if doPreplay
            figure
            tmp = (preVecSVMs-0.5);
            step = 0.01;
            cm = [[[0:step:1]'; ones(length(0:step:1),1)] ...
                [[0:step:1 1:-step:0]'] ...
                flipud([[0:step:1]'; ones(length(0:step:1),1)])];
            for gi = 1:length(envLabel)
                imagesc(nanmedian(tmp(:,:,ismember(envs,envLabel(gi))),3));
                colormap(cm);
                caxis([-0.5 0.5])
                saveFig(gcf,[root '/BurstVec_SVM_' envLabel{gi}],[{'tiff'} {'pdf'}])
            end
            
            doC = bsxfun(@times,ismember(envs,{'sq1'}),ismember(envs,{'g1'})');
            s1Tog1 = reshape(tmp(logical(repmat(doC, ...
                [1 1 length(sessions)]))),[nansum(doC(:)) length(sessions)]);
            doC = bsxfun(@times,ismember(envs,{'g1'}),ismember(envs,{'sq1'})');
            g1Tosq1 = -reshape(tmp(logical(repmat(doC, ...
                [1 1 length(sessions)]))),[nansum(doC(:)) length(sessions)]);
            comb = [s1Tog1; g1Tosq1];
            [a b] = nanmax(abs(comb),[],1);
            for i = 1:length(sessions)
                maxbias(i) = comb(b(i),i);
            end
            
            transitionPlot(crossSim,envs,doComps, ...
                [root '/ConditionComparison_PV']);
        end
        
        %%%%%%%%%%%%%%%%%% RDM ANALYSIS BEGIN %%%%%%%%%%%%%%
        
        features = [];
        for i = 1:length(envLabel)
            envI = find(ismember(envs,envLabel(i)));
            for k = 1:length(envI)
                features(envI(k),:) = [k i];
            end
        end
        lagMat = abs(bsxfun(@minus,[1:length(sim(:,1,1))],[1:length(sim(:,1,1))]'));
        iterMat = abs(bsxfun(@minus,features(:,1),features(:,1)'));
        envMat = abs(bsxfun(@minus,features(:,2),features(:,2)'));
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
        
        if ismember(mouse,{'AKCA130'})
            morphPoint = [4 4 5 5 5 5];
            vec = [];
            for gi = 1:6
                vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
            end
            attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
                [features(:,2)'>=vec(1:length(features(:,2)))]'));
        end
        
        if ismember(mouse,{'AKCA133'})
            morphPoint = [5 5 5 5 5 5];
            vec = [];
            for gi = 1:6
                vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
            end
            attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
                [features(:,2)'>=vec(1:length(features(:,2)))]'));
        end
        
        if ismember(mouse,{'AKCA135'})
            morphPoint = [5 5 5 5 5 5];
            vec = [];
            for gi = 1:6
                vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
            end
            attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
                [features(:,2)'>=vec(1:length(features(:,2)))]'));
        end
        
        RDMs = cat(3,lagMat,attractorMat,envMat);
        
        
%         tsim = nanmedian(sim,3);
%         t1s = tsim(1:round(length(tsim)./2),1:round(length(tsim)./2));
%         t1a = attractorMat(1:round(length(tsim)./2),1:round(length(tsim)./2));
%         t2s = tsim(round(length(tsim)./2)+1:end,round(length(tsim)./2)+1:end);
%         t2a = attractorMat(round(length(tsim)./2)+1:end,round(length(tsim)./2)+1:end);
%         firstHalfDecorr = nanmean(t1s(~logical(t1a)))-nanmean(t1s(logical(t1a)));
%         secondHalfDecorr = nanmean(t2s(~logical(t2a)))-nanmean(t2s(logical(t2a)));
        
%         fitAttractorDecorr(nanmedian(sim,3),RDMs(:,:,[1 2]));
        

%         classifyCellRDMs(sim,-RDMs,root);
        [r2penalty tr2] = classifyCellRDMswInteract(sim,-RDMs,root);
%         classifyCellRDMswInteract2D(sim,-RDMs(:,:,1:2),root);
%         [outPenalty tr2] = classifyCellRDMs2D(sim,-RDMs(:,:,1:2),root);
%         classifyCellRDMs(kDeformSim,-RDMs,[root '/GeoCorrected']);
%         classifyCellRDMswInteract(kDeformSim,-RDMs,[root '/GeoCorrected']);
%         classifyCellRDMswInteract2D(kDeformSim,-RDMs(:,:,1:2),[root '/GeoCorrected']);
%         classifyCellRDMs2D(kDeformSim,-RDMs(:,:,1:2),[root '/GeoCorrected']);


        valVsClass(nanmedian(ushc,2),r2penalty,tr2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}])
        saveFig(gcf,[root '/Classification/SHC_vs_tr2'],[{'pdf'} {'tiff'}])
        
        allR2Penalty = [allR2Penalty; r2penalty];
        allTR2 = [allTR2; tr2];
        
        
        allCellProps = [allCellProps; nanmedian(ushc,2) nanmedian(usic,2) nanmedian(umfr,2)];

%         if ismember(mouse,{'AKCA127'})
%             morphPoint = [4 4 4 4];
%             vec = [];
%             for gi = 1:4
%                 vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
%             end
%             attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
%                 [features(:,2)'>=vec(1:length(features(:,2)))]'));
%         end
        
%         lagMat = lagMat-nanmean(lagMat(:));
%         iterMat = iterMat-nanmean(iterMat(:));
%         envMat = envMat-nanmean(envMat(:));

%         sfpXVal(msfps,outPenalty,tr2./nanmax(tr2));
%         sfpXVal(msfps(:,:,tr2>0.1),outPenalty(tr2>0.1,:));

        
% % %         fits = nan(length(sim(1,1,:)),6);
% % %         subFits = nan(length(sim(1,1,:)),3);
% % %         for k = 1:length(sim(1,1,:))
% % %             tmp = sim(:,:,k);
% % %             if all(isnan(tmp(:)))
% % %                 continue
% % %             end
% % %             fits(k,:) = [corr(tmp(~isnan(tmp)),lagMat(~isnan(tmp)),'type','kendall') ...
% % %                 corr(tmp(~isnan(tmp)),attractorMat(~isnan(tmp)),'type','kendall') ...
% % %                 corr(tmp(~isnan(tmp)),envMat(~isnan(tmp)),'type','kendall') ...
% % %                 corr(tmp(~isnan(tmp)),attractorMat(~isnan(tmp)).*lagMat(~isnan(tmp)),'type','kendall') ...
% % %                 corr(tmp(~isnan(tmp)),envMat(~isnan(tmp)).*lagMat(~isnan(tmp)),'type','kendall') ...
% % %                 nansum(~isnan(tmp(:)))];
% % %             subTmp = tmp-fits(k,2).*attractorMat;
% % %             subFits(k,:) = [corr(subTmp(~isnan(tmp)),iterMat(~isnan(tmp)),'type','kendall') ...
% % %                 corr(subTmp(~isnan(tmp)),attractorMat(~isnan(tmp)),'type','kendall') ...
% % %                 corr(subTmp(~isnan(tmp)),envMat(~isnan(tmp)),'type','kendall')];
% % %             
% % %         end
% % %         fits = fits(:,[1 2 3 4 6]);
% % % 
% % %         minSamples = nchoosek(length(sessions),2).*0.25;
% % %         colors = bsxfun(@times,[0 0.25 1],-fits(fits(:,end)>minSamples,2));
% % %         colors = colors+bsxfun(@times,[1 0 0],-fits(fits(:,end)>minSamples,1));
% % %         colors(colors<0) = 0;
% % %         colors = colors.^(1./4);
% % %         colors = [colors./nanmax(colors(:))];
% % %         figure
% % %         set(gcf,'position',[50 50 350 350])
% % %         scatter(-fits(fits(:,end)>minSamples,2), ...
% % %             -fits(fits(:,end)>minSamples,1),...
% % %             30,colors,'filled')
% % %         hold on
% % %         xlabel('Attractor RDM Fit (Kendall''s \tau)')
% % %         ylabel('Time RDM Fit (Kendall''s \tau)')
% % %         set(gca,'xlim',[-0.4 0.8],'ylim',[-0.4 0.8])
% % %         axis square
% % %         plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
% % %         plot([0 0],get(gca,'ylim'),'linestyle','--','color','k')
% % %         saveFig(gcf,[root '/RDM_Analysis_Cellwise'],[{'pdf'} {'tiff'}])
% % %    
% % %         [goodInds] = find(fits(:,end)>minSamples);
% % %         figure
% % %         set(gcf,'position',[50 50 700 700])
% % %         subplot(2,2,1)
% % %         imagesc(-lagMat)
% % %         axis off
% % %         subplot(2,2,2)
% % %         imagesc(-attractorMat)
% % %         axis off
% % %         subplot(2,2,3)
% % %         [a best] = sort(abs(fits(fits(:,end)>minSamples,1)),'descend');
% % %         ex = squarify(sim(:,:,(goodInds(best(3)))));
% % %         imagesc(ex)
% % %         alpha(double(~isnan(ex)))
% % %         axis off
% % %         subplot(2,2,4)
% % %         [a best] = sort(abs(fits(fits(:,end)>minSamples,2)),'descend');
% % %         ex = squarify(sim(:,:,(goodInds(best(1)))));
% % %         imagesc(ex)
% % %         alpha(double(~isnan(ex)))
% % %         axis off
% % %         saveFig(gcf,[root '/RDM_Analysis_Cellwise_Examples'],[{'pdf'} {'tiff'}])
% % %         
% % %         %%%% Correlation with firing properties
% % %         %%% SHC vs Attractor
% % %         
% % %         %%%%%% Group by max within familiar environment
% % %         
% % % %         best_shc = nanmax([nanmedian(ushc(:,ismember(envs,{'sq1'})),2) ...
% % % %             nanmedian(ushc(:,ismember(envs,{'g1'})),2)],[],2);
% % %         valVsRDMFits(nanmean(ushc,2),fits,minSamples)
% % %         saveFig(gcf,[root '/RDM_Fits_vs_SHC'],[{'pdf'} {'tiff'}])
% % %         
% % % %         best_sic = nanmax([nanmedian(usic(:,ismember(envs,{'sq1'})),2) ...
% % % %             nanmedian(usic(:,ismember(envs,{'g1'})),2)],[],2);
% % %         valVsRDMFits(nanmean(usic,2),fits,minSamples)
% % %         saveFig(gcf,[root '/RDM_Fits_vs_SIC'],[{'pdf'} {'tiff'}])
% % % 
% % %         valVsRDMFits(nanmean(upfs,2),fits,minSamples)
% % %         saveFig(gcf,[root '/RDM_Fits_vs_PlaceFieldSize'],[{'pdf'} {'tiff'}])
% % %         
% % %         valVsRDMFits(nanmean(umfr,2),fits,minSamples)
% % %         saveFig(gcf,[root '/RDM_Fits_vs_MFRs'],[{'pdf'} {'tiff'}])
% % % %         
% % % %         
% % % %         plotStackConMaps(permute(um(:,:,goodInds(best(3)),:),[1 2 4 3]),envs,envLabel)
%         
        %%%%%%%%%%%%%%%%% RDM ANALYSIS END %%%%%%%%%%%%%%%%%
        
        %%%%%%%%% GEO FITS VERSUS DRIFT / ATTRACTOR LOADINGS
        
%         figure
%         set(gcf,'position',[50 50 1400 700])
%         for i =1:4
%             geofits{i} = permute(nansum(nansum(kDeformFit==i,1),2)./ ...
%                 nansum(nansum(~isnan(kDeformFit),1),2),[3 4 1 2]);
%             subplot(2,4,i)
%             scatter(geofits{i}(fits(:,end)>minSamples),-fits(fits(:,end)>minSamples,1))
%             lsline
%             subplot(2,4,i+4)
%             scatter(geofits{i}(fits(:,end)>minSamples),-fits(fits(:,end)>minSamples,2))
%             lsline
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         mds3D(sim,envs,envLabel,[root '/MDS_3D.gif']);
%         mds3D(sim(:,:,alignedSessions<=16),envs,envLabel,[root '/MDS_3D_LessRegistered.gif']);
%         mds3D(sim(:,:,alignedSessions>16),envs,envLabel,[root '/MDS_3D_MoreRegistered.gif']);
%         
%         mkSessionComps(sim,envs,doComps,[root '/ConditionComparison_Cellwise']);
%         step = 5;
%         for i = step:step:35
%             mkSessionComps(sim(:,:,(alignedSessions>(i-step))&(alignedSessions<=(i))),envs,doComps, ...
%                 [root '/ConditionComparison_Cellwise_Registered_' num2str(i-step) '_to_' num2str(i)]);
%         end
%         tmp = repmat({[]},[size(kDeformSim(:,:,1))]);
%         for i = 1:length(kDeformSim(1,:,1))
%             for j = 1:length(kDeformSim(1,:,1))
%                 tmp{i,j} = permute(sim(i,j,~isnan(kDeformSim(i,j,:))),[3 1 2]);
%             end
%         end
% 
%         for gi = 1:6:length(kDeformSim(1,:,1))
%             params = transitionPlot(tmp(gi:nanmin(gi+7,length(kDeformSim(1,:,1))), ...
%                 gi:nanmin(gi+7,length(kDeformSim(1,:,1)))), ...
%                 envs(gi:nanmin(gi+7,length(kDeformSim(1,:,1)))),doComps, ...
%                 [root '/ConditionComparison_Correlations_Cellwise_Sequence' num2str(((gi-1)./6) + 1)]);
%         end
    end
    
    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
    for ti = 0:0.1:0.5
        plotWeightClass(allR2Penalty(~isnan(allTR2)&allTR2>ti,:),allTR2(allTR2>ti),...
            [root '/OVERALL/Dropout_withInteractions_r2_threshold_' num2str(ti) '.gif'],false);
    end

    valVsClass(allCellProps(:,1),allR2Penalty,allTR2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}])
    saveFig(gcf,[root '/OVERALL/shc_vs_tr2'],[{'pdf'} {'tiff'}])
    
%     valVsClass(allCellProps(:,3),allR2Penalty,allTR2, ...
%             [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}])
%     saveFig(gcf,[root '/OVERALL/mfr_vs_tr2'],[{'pdf'} {'tiff'}])
end




















