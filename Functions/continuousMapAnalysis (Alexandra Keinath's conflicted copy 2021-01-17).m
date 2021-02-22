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
    allClassXTime = []; %repmat({[]},[4 1]);
    allClassXSeq = [];
    compAngDiffs = [];
    sigParams = [];
    classRSMs = repmat({[]},[2 3]);
    allDecayParams = repmat({[]},[1 3]);
    allDecayFunctions = repmat({[]},[1 3]);
    allClassDurs = repmat({[]},[1 4]);
    allClassFirstLast = repmat({[]},[1 4]);
    
    shuffle_r2penalty = [];
    shuffle_tr2 = [];
    
    envSize = [17 17];
    includeLags = inf;
    pThresh = 0.01; %0.025;
    tr2Thresh = 0.1915;
    
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
        close all
        drawnow;
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)
            if doAnatomy
                s = load(sessions{si},'processed','exclude');
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
%             if doAnatomy
%                 MeanFrames{si} = s.calcium.meanFrame;
%             end

            
%             figure(1)
%             set(gcf,'position',[50 50 1600 600])
%             subplot(4,8,si)
%             imagesc(s.processed.meanFrame);
%             axis equal
%             axis off
%             colormap gray
%             
%             figure(2)
%             set(gcf,'position',[50 50 1600 600])
%             subplot(4,8,si)
%             tmp = s.processed.SFPs(:,:,s.processed.exclude.SFPs);
%             tmp = tmp./repmat(nanmax(nanmax(tmp,[],1),[],2),size(tmp(:,:,1)));
%             imagesc(nanmax(tmp,[],3).^3);
%             axis equal
%             axis off
        end  
%         figure(1)
%         saveFig(gcf,[root '/MeanFrames'],[{'pdf'} {'tiff'}])
%         figure(2)
%         saveFig(gcf,[root '/SFPs'],[{'pdf'} {'tiff'}])
%         close all
%         drawnow;
%         durat = toc;
%         fprintf([num2str(durat) ' s\n']);
        
        minVecTimes = nanmin(cellfun(@size,simVecs,repmat({2},size(simVecs))));
        uvecSim = nan([length(alignMap{1}(:,1)) minVecTimes length(sessions)]);
        umfr = nan([length(alignMap{1}(:,1)) length(sessions)]);
        upfs = nan([length(alignMap{1}(:,1)) length(sessions)]);
        usic = nan([length(alignMap{1}(:,1)) length(sessions)]);
        ushc = nan([length(alignMap{1}(:,1)) length(sessions)]);
        uclass = nan([length(alignMap{1}(:,1)) length(sessions)]);
        uIsPC = false([length(alignMap{1}(:,1)) length(sessions)]);
        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        
        for si = 1:length(sessions)
            
            umfr(alignMap{1}(:,si)~=0,si) = MFRs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            upfs(alignMap{1}(:,si)~=0,si) = PFSs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            usic(alignMap{1}(:,si)~=0,si) = SICs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            ushc(alignMap{1}(:,si)~=0,si) = SHCs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            uIsPC(alignMap{1}(:,si)~=0,si) = isPC{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
            
        end
        
        
        if doAnatomy
            for si = 1:length(sessions)
                SFPs{si} = permute(SFPs{si},[2 3 1]);
            end
            fprintf('\t\tCompressing Spatial Footprints... ');
            tic
            usfps = nan([size(SFPs{1}(:,:,1)) length(sessions)]);
            for k = 1:length(alignMap{1}(:,1))
                tmp = nan([size(SFPs{1}(:,:,1)) length(sessions)]);
                for si = 1:length(sessions)
                    if alignMap{1}(k,si)~=0
                        tmp(:,:,si) = SFPs{si}(:,:,alignMap{1}(k,si));
                    end
                end
                usfps(:,:,k) = nanmean(tmp,3);
                usfps(:,:,k) = usfps(:,:,k)./nanmax(nanmax(usfps(:,:,k)));
            end
            durat = toc;
            fprintf([num2str(durat) ' s\n']);
        end
        
        if doAnatomy
            msfps = nanmean(usfps,4);
        end
      
        
        sim = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
        kCrossCOM = repmat({[]},length(sessions));
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing pairwise map comparisons: '])
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
                fprintf(str);
                strLength = length(str);
                
                tmp1 = um(:,:,:,si);
                tmp2 = um(:,:,:,sj);
            
                sim(si,sj,:) = xcorr3transform(tmp1,tmp2,[0 0 0]);

            end
        end       
        
%         %%%% shuffle rate maps
%         Sum = um;
%         for k = 1:length(um(1,1,:,1))
%             goodI = find(permute(~all(all(isnan(um(:,:,k,:)),1),2),[4 1 2 3]));
%             Sum(:,:,k,goodI) = um(:,:,k,goodI(randperm(length(goodI))));
%         end
%         
%         ssim = nan([length(sessions) length(sessions) length(alignMap{1}(:,1))]);
%         iter = 0;
%         strLength = 0;
%         fprintf(['\n\t\tShuffled within cell: '])
%         for si = 1:length(sessions)
%             for sj = si+1:length(sessions)
%                 
%                 iter = iter+1;
%                 fprintf(repmat('\b',[1 strLength]));
%                 str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
%                 fprintf(str);
%                 strLength = length(str);
%                 
%                 tmp1 = Sum(:,:,:,si);
%                 tmp2 = Sum(:,:,:,sj);
%             
%                 ssim(si,sj,:) = xcorr3transform(tmp1,tmp2,[0 0 0]);
% 
%             end
%         end
        
        if doCOM
            [a condGroup] = ismember(envs,envLabel);
            comVecPlot(groupMat(kCrossCOM,condGroup,false),[root '/COM_Vectors']);
        end


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

        if doMDS
            compAngDiffs = [compAngDiffs; mds2D(sim,envs,envLabel,[root '/MDS_2D'])];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEQUENCE ANALYSIS
        
        morphPoint = [];
        blah = find(ismember(upiece{mi},'/'),2,'last');
        mouseID = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end);
        for i = 1:floor(length(envs)./6)
            goodS = ((i-1).*6)+1:nanmin(((i-1).*6)+8,length(envs));
            tv = sim(goodS,goodS,:);
            tpc = uIsPC(:,goodS);
            ttv = repmat({[]},[size(tv(:,:,1))]);
            for j = 1:length(ttv)
                for k = 1:length(ttv)
%                     ttv{j,k} = permute(tv(j,k,:),[3 1 2]);
                    
                    ttv{j,k} = permute(tv(j,k,tpc(:,j)|tpc(:,k)),[3 1 2]);
                end
            end
            
            params = transitionPlot(ttv,envs(goodS),doComps, ...
                [root '/SequenceAnalyses/AllPlaceCells/Sequence_' num2str(i)]);
            
            if ~isempty(params)
                sigParams = [sigParams; params.sigmoidal_mad ...
                    params.sigmoidal_intercept-3.5 params.sigmoidal_max_sep mi i];
                morphPoint = [morphPoint ceil(params.sigmoidal_intercept)];
            end
        end
        morphPoint(end+1) = morphPoint(end);
        
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
        
        vec = [];
        for gi = 1:6
            vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
        end
        attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
            [features(:,2)'>=vec(1:length(features(:,2)))]'));
        
        mouse = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end);

        
        RDMs = cat(3,lagMat,attractorMat,envMat);
      
        %%%% make rdm example figure
        
        RDMs = RDMs./repmat(nanmax(nanmax(RDMs,[],1),[],2),[size(RDMs(:,:,1))]);
        RDMs = 1-RDMs;
        
%         figure
%         set(gcf,'position',[50 50 900 300])
%         subplot(1,3,1)
%         image(cat(3,RDMs(:,:,1),zeros(size(RDMs(:,:,1))),zeros(size(RDMs(:,:,1)))))
% %         axis equal
%         axis square
%         axis off
%         subplot(1,3,2)
%         image(cat(3,zeros(size(RDMs(:,:,1))),zeros(size(RDMs(:,:,1))),RDMs(:,:,2)))
% %         axis equal
%         axis square
%         axis off
%         subplot(1,3,3)
%         image(cat(3,zeros(size(RDMs(:,:,1))),RDMs(:,:,3),zeros(size(RDMs(:,:,1)))))
% %         axis equal
%         axis square
%         axis off
%         
%         [a b] = sort(tr2,'descend');
%         figure
%         set(gcf,'position',[50 50 300 300])
%         tmp = squarify(sim(:,:,b(375)));
%         imagesc(tmp)
%         alpha(double(~isnan(tmp)))
%         axis square
%         axis off

        

        [r2penalty tr2] = classifyCellRDMswInteract(sim,-RDMs,root);
        plotWeightClass(r2penalty(tr2>tr2Thresh,:),tr2(tr2>tr2Thresh), ...
            [root '/Classification/Dropout_withInteractions_r2.gif'],false);
        
        valVsClass(nanmedian(ushc,2),r2penalty,tr2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}]);
        saveFig(gcf,[root '/Classification/SHC_vs_tr2'],[{'pdf'} {'tiff'}]);
        
        allR2Penalty = [allR2Penalty; r2penalty];
        allTR2 = [allTR2; tr2];
        
        
        [ta tb] = classifyCellRDMswInteract(sim,-RDMs,root,true);
        
%         [ta tb] = classifyCellRDMswInteract(ssim,-RDMs,root);
        
        shuffle_r2penalty = [shuffle_r2penalty; ta];
        shuffle_tr2 = [shuffle_tr2; tb];
        
        allCellProps = [allCellProps; nanmedian(ushc,2) nanmedian(usic,2) nanmedian(umfr,2)];

        [a b] = nanmax(r2penalty,[],2);
        b(all(isnan(r2penalty),2)) = nan;
        cellClass = b;
        cellClass(tr2<tr2Thresh) = 0;
        
        %%%%% Fit decay of time-loading cells
        for i = 1:3
            [ap val] = fitTimeDecay(sim(:,:,b==i&tr2>tr2Thresh));
            allDecayParams{i} = [allDecayParams{i}; ap];
            allDecayFunctions{i} = [allDecayFunctions{i}; val];
        end
        
        uclass = bsxfun(@times,~isnan(umfr),b);
        uclass(tr2<tr2Thresh,:) = 0;
        uclass(isnan(umfr)) = nan;
        
        tmp = nan(4,length(uclass(1,:)));
        
        classDurs = repmat({[]},[1 4]);
        for i = 1:4
            goodClass = any(uclass==i-1,2);
            classDurs{i} = nansum(~isnan(uclass(goodClass,:)),2);
            allClassDurs{i} = [allClassDurs{i};  classDurs{i}];
            
            fl = nan(nansum(goodClass),2);
            igc = find(goodClass);
            for k = igc'
                blah = ~isnan(uclass(k,:));
                fl(igc==k,:) = [find(blah,1,'first') find(blah,1,'last')];
            end
            allClassFirstLast{i} = [allClassFirstLast{i};  fl ]; %fl(:,2)-fl(:,1)];
            for j = 1:length(uclass(1,:))
                tmp(i,j) = nansum(uclass(:,j)==(i-1))./nansum(~isnan(uclass(:,j)));
            end
        end
        
        allClassXTime = [allClassXTime {tmp}];
        
        classXSeq = [];
        for i = 1:floor(length(envs)./6)
            goodS = ((i-1).*6)+1:nanmin(((i-1).*6)+8,length(envs));
            classXSeq = [classXSeq; nanmedian(tmp(:,goodS),2)'];
        end
        allClassXSeq = [allClassXSeq; classXSeq];
        
        
        %%%%%%%%%%%JUST CELL GROUPS SEQUENCE ANALYSIS
        %%%%%%%%%%%%%%%%%%%%% Subpopulation RDMs %%%%%%%%%%%%%%%%%%%%%%
        
        labels = [{'TimeCells'} {'AttractorCells'} {'GeometryCells'}];
        
        lims = [1 -1];
        for group = 1:3
            tmp = nanmedian(sim(:,:,cellClass==group),3);
            lims = [nanmin(lims(1),nanmin(tmp(:))) ...
                nanmax(lims(2),nanmax(tmp(:)))];
        end
        lims = [floor(lims(1).*10)./10 ceil(lims(2).*10)./10];
        
        figure
        set(gcf,'position',[50 50 900 600]);
        for group = 1:3
            
            classRSMs{1,group} = cat(3,classRSMs{1,group},sim(:,:,cellClass==group));
            
            subplot(2,3,group)
            tmp = squarify(nanmedian(sim(:,:,cellClass==group),3));
            imagesc(tmp)
            alpha(double(~isnan(tmp)))
            caxis(lims)
            axis square
            axis equal
            axis off
            colorbar
            
            newOrder = [];
            for i = 1:length(envLabel)
                newOrder = [newOrder; find(ismember(envs,envLabel(i)))];
            end
            tmp = tmp(newOrder,newOrder);    
            
            classRSMs{2,group} = cat(3,classRSMs{2,group},tmp);
            
            subplot(2,3,group+3)
            imagesc(tmp)
            alpha(double(~isnan(tmp)))
            caxis(lims)
            axis square
            axis equal
            axis off
            colorbar
        end
        saveFig(gcf,[root '/RDMs_Subpopulations'],[{'pdf'} {'tiff'}])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%         for group = 1:3
%             
%             for i = 1:floor(length(envs)./6)
%                 goodS = ((i-1).*6)+1:nanmin(((i-1).*6)+8,length(envs));
%                 tv = sim(goodS,goodS,cellClass==group);
%                 tpc = uIsPC(cellClass==group,goodS);
%                 ttv = repmat({[]},[size(tv(:,:,1))]);
%                 tpc = true(size(tpc));
%                 for j = 1:length(ttv)
%                     for k = 1:length(ttv)
%                         ttv{j,k} = permute(tv(j,k,tpc(:,j)|tpc(:,k)),[3 1 2]);
%                     end
%                 end
% 
%                 params = transitionPlot(ttv,envs(goodS),doComps, ...
%                     [root '/SequenceAnalyses/' labels{group} '/Sequence_' num2str(i)]);
%             end
%         end
        
        close all
        drawnow
%         
%         %%%%%%%% Plot maps with those with high tr2
%         for k = 1:length(tr2)
%             if tr2(k) > tr2Thresh
%                 plotStackConMaps(um(:,:,k,:),envs);
%                 saveFig(gcf,[root '/ContinuousMaps/TimeOrdered/' labels{cellClass(k)} '/Cell_' num2str(k)],[{'pdf'} {'tiff'}]);
%                 plotStackConMaps(um(:,:,k,:),envs,envLabel);
%                 saveFig(gcf,[root '/ContinuousMaps/GeometryOrdered/' labels{cellClass(k)} '/Cell_' num2str(k)],[{'pdf'} {'tiff'}]);
%                  close all
%                 drawnow
%             end
%         end
        
        if doAnatomy
            toPlot = cat(3,nanmax(usfps(:,:,cellClass==1).^8,[],3), ...
                nanmax(usfps(:,:,cellClass==3).^8,[],3), ...
                nanmax(usfps(:,:,cellClass==2).^8,[],3)).*255;
            image(toPlot)
            saveFig(gcf,[root '/SFP_x_class'],[{'pdf'} {'tiff'}])
        end
        
        fprintf('\n');
    end
    
%     tmp = shuffle_tr2;
%     tmp = sort(tmp,'descend');
%     tmp(isnan(tmp)) = [];
%     fid = fopen('Threshold_tr2.txt','w');
%     for i = [0.01 0.05 0.1 0.2 0.3]
%         val = tmp(round(length(tmp).*i));
%         fprintf(fid,'\n%.0f Percentile: %0.4f',[i.*100 val]);
%     end
%     fclose(fid);
% 
%     figure
%     set(gcf,'position',[50 50 250 250])
%     cumHist([{allTR2} {shuffle_tr2}],[0:0.01:1])
%     xlabel('Full model tr2')
%     ylabel('Cumulative proportion')
%     saveFig(gcf,[root '/Overall/tr2_vs_shuffled'],[{'pdf'} {'tiff'}])

    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
    
    figure
    set(gcf,'position',[50 50 900 600]);
    for group = 1:3
        subplot(2,3,group)
        tmp = squarify(nanmedian(classRSMs{1,group},3));
        imagesc(tmp)
        alpha(double(~isnan(tmp)))
        axis square
        axis equal
        axis off
        caxis([-0.2 0.9])
        colorbar


        tmp = squarify(nanmedian(classRSMs{2,group},3));
        
        subplot(2,3,group+3)
        imagesc(tmp)
        alpha(double(~isnan(tmp)))
        axis square
        axis equal
        axis off
        caxis([-0.2 0.9])
        colorbar
    end
    saveFig(gcf,[root '/Overall/RDMs_Subpopulations' labels{group}],[{'pdf'} {'tiff'}])
    
    %%%%%%%%%%%%%%%%%% Temporal decay analysis
    toPlot = [{[]} {[]} {[]}];
    for i = 1:3
        toPlot{i} = allDecayParams{i}(:,2); %./allDecayParams{i}(:,1);
    end
    a = cellfun(@nanmean,toPlot);
    b = cellfun(@nanstd,toPlot)./sqrt(cellfun(@numel,toPlot));
    c = [a; b];
    fid = fopen('Stats_DecayAnalysis.txt','w');
    fprintf(fid,'\n\tM+/-SEM: %0.4f +/- %0.6f',c);
    fprintf(fid,'\n\n\t\tSigned Rank vs 0\n');
    for i = 1:3
        [pval h stats] = signrank(toPlot{i},0);
        fprintf(fid,['\n(Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
            [i stats.signedrank stats.zval pval]);
    end
    fprintf(fid,'\n\n\t\tRank sum pairwise comparisons\n');
    for i = 1:3
        for j = i+1:3
            [pval h stats] = ranksum(toPlot{i},toPlot{j});
            fprintf(fid,['\n(Group %.0f vs Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
                [i j stats.ranksum stats.zval pval]);
        end
    end
    fclose(fid);
    
    figure
    set(gcf,'position',[50 50 225 250],'color','w')
    mkBow(toPlot,[{'Time'} {'Attractor'} {'Geometry'}]);
    hold on
    plot(get(gca,'xlim'),[0 0],'linestyle','--')
    set(gca,'outerposition',[0 0 1 0.9])
    ylabel('Similarity x lag (r/day)')
    xlabel('Cell group')
    saveFig(gcf,[root '/Overall/Decay_Slope'],[{'pdf'} {'tiff'}])
    
%     figure
%     mkLine(allDecayFunctions,[0:1:31])
    
%     figure
%     set(gcf,'position',[50 50 900 300])
%     subplot(1,3,1)
%     mkBow(allDecayFunctions{1}(:,2:end))
%     subplot(1,3,2)
%     mkBow(allDecayFunctions{2}(:,2:end))
%     subplot(1,3,3)
%     mkBow(allDecayFunctions{2}(:,2:end))
    
    %%%%%%%%%%%%%%%%%%%%% SEQUENCE ANALYSES %%%%%%%%%%%%%%%%%%%%%
    
    figure
    set(gcf,'position',[50 50 600 300])
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,3);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    subplot(1,2,1)
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
%     plot(nanmean(toPlot)','marker','none','color','k','linewidth',2,'linestyle','-')
    set(gca,'xlim',[0 6])
%     h = mkGraph(toPlot);
    ylabel('Maximum decorrelation (Z(\Deltar))')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'ylim',[0 1.2])
    subplot(1,2,2)
    tmp = bsxfun(@minus,toPlot(:,2:end),toPlot(:,1));
    h = plot([zeros(length(toPlot(:,1)),1) tmp]', ...
        'marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
%     mkGraph(bsxfun(@minus,toPlot(:,2:end),toPlot(:,1)))
    set(gca,'xlim',[0 6])
%     h = mkGraph(toPlot);
    ylabel('Maximum decorrelation (Z(\Deltar))')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'ylim',[-0.4 0.4])
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    saveFig(gcf,[root ...
        '/Overall/SequenceAnalyses/Maximum_Separation'],[{'pdf'} {'tiff'}]);
    
    pvals = [];
    zstats = [];
    for i = 1:4
        [pval h stats] = signrank(tmp(:,i),0);
        pvals = [pvals pval];
        zstats = [zstats stats.signedrank];
    end
    
%     pvals = [];
%     zstats = [];
%     for i = 1:4
%         [h pval ci stats] = ttest(tmp(:,i),0);
%         pvals = [pvals pval];
%         zstats = [zstats stats.tstat];
%     end
    
    
    figure
    set(gcf,'position',[50 50 300 300])
    hold on
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,2);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
    set(gca,'xlim',[0 6])
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    ylabel('median(|r_{sq1} - r_{g1}|)')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'ylim',[-2.5 2.5],'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('Transition Point')
    xlabel('Morph Sequence Number')
    saveFig(gcf,[root ...
        '/Overall/SequenceAnalyses/SigFitIntercept'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 50 600 300])
    hold on
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,1);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    subplot(1,2,1)
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
    plot(nanmean(toPlot)','marker','none','color','k','linewidth',2,'linestyle','-')
    set(gca,'xlim',[0 6])
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    ylabel('median(|r_{sq1} - r_{g1}|)')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('mean(abs(sig-lin))')
    xlabel('Morph Sequence Number')
    subplot(1,2,2)
    h = plot([zeros(length(toPlot(:,1)),1) bsxfun(@minus,toPlot(:,2:end),toPlot(:,1))]', ...
        'marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    set(gca,'xlim',[0 6])
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    ylabel('median(|r_{sq1} - r_{g1}|)')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('mean(abs(sig-lin))')
    xlabel('Morph Sequence Number')
    
    saveFig(gcf,[root ...
        '/Overall/SequenceAnalyses/MeanAbsDiffBtwn_SigLin'],[{'pdf'} {'tiff'}]);
    
    
    labels = [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}];
    figure
    set(gcf,'position',[50 50 900 300])
    for i = 1:3
        subplot(1,3,i)
        scatter(allClassXSeq(:,i+1),sigParams(:,3))
        xlabel(['Proportion ' lower(labels{i})]);
        ylabel('Maximum Map Separation');
        set(gca,'xlim',[0 0.4],'ylim',[0 1])
        lsline
        
        [rval pval] = corr(allClassXSeq(:,i+1),sigParams(:,3));
        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.9;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.1;
        text(textX,textY,sprintf(['r = %.3f\np = %.2e'],[rval pval]), ...
            'fontweight','normal','fontname','arial','fontsize',10);
    end
    saveFig(gcf,[root ...
        '/Overall/SequenceAnalyses/Maximum_Separation_vs_ProportionClass'],[{'pdf'} {'tiff'}]);
    
%     figure
%     set(gcf,'position',[50 50 1200 300])
%     for i = 1:3
%         subplot(1,4,i)
%         scatter(allClassXSeq(:,i+1),sigParams(:,1))
%         xlabel(['Proportion ' lower(labels{i})]);
%         ylabel('mean(|sig-lin|)');
%         lsline
%         
%         [rval pval] = corr(allClassXSeq(:,i+1),sigParams(:,1));
%         textX = get(gca,'xlim');
%         textX = textX(2) - [textX(2)-textX(1)].*0.9;
%         textY = get(gca,'ylim');
%         textY = textY(2) - [textY(2)-textY(1)].*0.1;
%         text(textX,textY,sprintf(['r = %.3f\np = %.2e'],[rval pval]), ...
%             'fontweight','normal','fontname','arial','fontsize',10);
%     end
%     
%     subplot(1,4,4)
%     scatter(allClassXSeq(:,3)-allClassXSeq(:,4),sigParams(:,1))
%     xlabel(['Proportion attractor - geometry cells']);
%     ylabel('mean(|sig-lin|)');
%     lsline
% 
%     [rval pval] = corr(allClassXSeq(:,3)-allClassXSeq(:,4),sigParams(:,1));
%     textX = get(gca,'xlim');
%     textX = textX(2) - [textX(2)-textX(1)].*0.9;
%     textY = get(gca,'ylim');
%     textY = textY(2) - [textY(2)-textY(1)].*0.1;
%     text(textX,textY,sprintf(['r = %.3f\np = %.2e'],[rval pval]), ...
%         'fontweight','normal','fontname','arial','fontsize',10);
%     
%     saveFig(gcf,[root ...
%         '/Overall/SequenceAnalyses/MeanAbsDiff_vs_ProportionClass'],[{'pdf'} {'tiff'}]);
    
    
    %%%%%%%%%%% Ranked factor plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ti = [0:0.05:0.3 tr2Thresh]
    
        
        allStats = [];
        ta = sort(allR2Penalty,2,'descend');
        sa = sort(shuffle_r2penalty,2,'descend');
        toPlot = repmat({[]},[2 3]);
        for i = 1:3
            toPlot{1,i} = ta(~isnan(allTR2)&allTR2>ti,i);
            toPlot{2,i} = sa(~isnan(shuffle_tr2)&shuffle_tr2>ti,i);
            
            a = toPlot{1,i}(~isnan(toPlot{1,i}));
            b = toPlot{2,i}(~isnan(toPlot{2,i}));
            [pval h stats] = ranksum(a,b);
            allStats = [allStats [pval stats.zval]'];
        end
        
        figure
        set(gcf,'position',[50 50 225 250],'color','w')
        mkBow(toPlot,[{'1st'} {'2nd'} {'3rd'}]);
        xlabel('Factor rank')
        ylabel('Attributable r^2')
        for i = 1:3
            text(i,1.11,sprintf('Z=%.2f\np=%.1e\n', ...
                abs(allStats(2,i)),allStats(1,i)),...
                'horizontalalignment','center','color','k','fontsize',7);
        end
        text(0.1,1.15,sprintf('Threshold:\nr2 > %0.2f',ti),...
            'horizontalalignment','center','color','k','fontsize',7);
        set(gca,'outerposition',[0 0 1 0.9])
        saveFig(gcf,[root '/Overall/Ranked_factor_tr2_threshold_' num2str(ti)],[{'pdf'} {'tiff'}])
        
        fid = fopen(['Stats_RankedFactor.txt'],'w');
        for i = 1:3
            [pval h stats] = ranksum(toPlot{1,i},toPlot{2,i});
            fprintf(fid,['\n(Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
                [i stats.ranksum stats.zval pval]);
        end
        fclose(fid);
    end
    
%     for ti = [0 tr2Thresh]
%         allStats = [];
%         [ta torder] = sort(allR2Penalty,2,'descend');
%         [sa sorder] = sort(shuffle_r2penalty,2,'descend');
%         toPlot = repmat({[]},[2 9]);
%         for i = 1:3
%             for j = 1:3
%                 toPlot{1,(j-1).*3+i} = allR2Penalty(~isnan(allTR2) & ...
%                     allTR2>ti & torder(:,1)==i,j);
%                 toPlot{2,(j-1).*3+i} = ...
%                     shuffle_r2penalty(~isnan(shuffle_tr2) & ...
%                     shuffle_tr2>ti&sorder(:,1)==i,j)
% 
%                 a = toPlot{1,i}(~isnan(toPlot{1,i}));
%                 b = toPlot{2,i}(~isnan(toPlot{2,i}));
%                 [pval h stats] = ranksum(a,b);
%                 allStats = [allStats [pval stats.ranksum]'];
%             end
%         end
%         
%         figure
%         set(gcf,'position',[50 50 225 300],'color','w')
%         mkBow(toPlot);
% %         xlabel('Factor rank')
% %         ylabel('Attributable r^2')
% %         for i = 1:3
% %             text(i,1.12,sprintf('Z=%.1e\np=%.1e',abs(allStats(2,i)),allStats(1,i)),...
% %                 'horizontalalignment','center','color','k','fontsize',7)
% %         end
% %         set(gca,'outerposition',[0 0 1 0.9])
% %         saveFig(gcf,[root '/Overall/Ranked_factor_tr2_threshold_' num2str(ti)],[{'pdf'} {'tiff'}])
%     end
    
    %%%%%%%%%% Duration characteristics x class

    figure
    set(gcf,'position',[50 50 225 250],'color','w')
    mkBow(allClassDurs(2:4),[{'Time'} {'Attractor'} {'Geometry'}]);
    set(gca,'outerposition',[0 0 1 0.9],'ylim',[0 32])
    ylabel('Number of tracked sessions')
    xlabel('Cell group')
    saveFig(gcf,[root '/Overall/Class_x_NumTrackedSessions'],[{'pdf'} {'tiff'}])
    
    toPlot = allClassDurs(2:4);
    fid = fopen('Stats_TrackedSessions.txt','w');
    for i = 1:3
        [pval h stats] = signrank(toPlot{i},0);
        fprintf(fid,['\n(Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
            [i stats.signedrank stats.zval pval]);
    end
    fprintf(fid,'\n\n\t\tRank sum pairwise comparisons\n');
    for i = 1:3
        for j = i+1:3
            [pval h stats] = ranksum(toPlot{i},toPlot{j});
            fprintf(fid,['\n(Group %.0f vs Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
                [i j stats.ranksum stats.zval pval]);
        end
    end
    fclose(fid);
    
%     tmpA = [{allClassFirstLast{1}(:,1)} {allClassFirstLast{2}(:,1)} ...
%         {allClassFirstLast{3}(:,1)} {allClassFirstLast{4}(:,1)}];
%     tmpB = [{allClassFirstLast{1}(:,2)} {allClassFirstLast{2}(:,2)} ...
%         {allClassFirstLast{3}(:,2)} {allClassFirstLast{4}(:,2)}];
%     
%     figure
%     set(gcf,'position',[50 50 400 400])
%     mkBowSessions(cellfun(@diff,allClassFirstLast, ...
%         [{[]} {[]} {[]} {[]}],[{[2]} {[2]} {[2]} {[2]}],'uniformoutput',false))
%     saveFig(gcf,[root '/Overall/Class_x_FLDiff'],[{'pdf'} {'tiff'}])
%     
%     figure
%     set(gcf,'position',[50 50 400 400])
%     mkBowSessions(tmpA)
%     saveFig(gcf,[root '/Overall/Class_xFirstSession'],[{'pdf'} {'tiff'}])
%     
%     figure
%     set(gcf,'position',[50 50 400 400])
%     mkBowSessions(tmpB)
%     saveFig(gcf,[root '/Overall/Class_xLastSession'],[{'pdf'} {'tiff'}])
    
    %%%%%%%%%%%%%%%%%%%%%
    
    for ti = [0 tr2Thresh]
        plotWeightClass(allR2Penalty(~isnan(allTR2)&allTR2>ti,:),allTR2(allTR2>ti),...
            [root '/Overall/Dropout_withInteractions_r2_threshold_' num2str(ti) '.gif'],false);
    end

    toPlot = repmat({[]},[1 3]);
    [a b] = nanmax(allR2Penalty,[],2);
    b(isnan(allTR2)) = nan;
    for i = 1:3
        toPlot{i} = allCellProps(b==i & allTR2 > tr2Thresh,1);
    end
    cellfun(@nanmean,toPlot)
    cellfun(@nanstd,toPlot)./sqrt(cellfun(@numel,toPlot))
    figure
    set(gcf,'position',[50 50 225 250],'color','w')
    mkBow(toPlot,[{'Time'} {'Attractor'} {'Geometry'}]);
    set(gca,'outerposition',[0 0 1 0.9])
    xlabel('Cell group')
    ylabel('Median within-session SHC (r)')
    saveFig(gcf,[root '/Overall/WithinSessionSHC'],[{'pdf'} {'tiff'}])
    
    a = cellfun(@nanmean,toPlot);
    b = cellfun(@nanstd,toPlot)./sqrt(cellfun(@numel,toPlot));
    c = [a; b];
    fid = fopen('Stats_WithinSessionSHC.txt','w');
    fprintf(fid,'\n\tM+/-SEM: %0.4f +/- %0.6f',c);
    fprintf(fid,'\n\n\t\tSigned Rank vs 0\n');
    for i = 1:3
        [pval h stats] = signrank(toPlot{i},0);
        fprintf(fid,['\n(Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
            [i stats.signedrank stats.zval pval]);
    end
    fprintf(fid,'\n\n\t\tRank sum pairwise comparisons\n');
    for i = 1:3
        for j = i+1:3
            [pval h stats] = ranksum(toPlot{i},toPlot{j});
            fprintf(fid,['\n(Group %.0f vs Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
                [i j stats.ranksum stats.zval pval]);
        end
    end
    fclose(fid);
    
    
    
    allStats = [];
    for i = 1:3
        for j = i+1:3
            a = toPlot{1,i};
            b = toPlot{1,j};
            [pval h stats] = ranksum(a,b);
            allStats = [allStats [pval stats.ranksum]'];
        end
    end
    
    
    atab = valVsClass(allCellProps(:,1),allR2Penalty,allTR2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}]);
    saveFig(gcf,[root '/Overall/shc_vs_tr2'],[{'pdf'} {'tiff'}]);
    cell2file(atab,[root '/Overall/shc_vs_tr2.txt']);
    
    valVsClass(allCellProps(:,3),allR2Penalty,allTR2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}]);
    saveFig(gcf,[root '/Overall/mfr_vs_tr2'],[{'pdf'} {'tiff'}]);

    fid = fopen('Stats_Embedded Angles.txt','w');
    fprintf(fid,'\n\n\t**************************************\n\n');
    fprintf(fid,'\tEmbedded Angle absolute differences:');
    fprintf(fid,'\n\t\t%0.2f',abs(compAngDiffs));
    [pval zval] = circ_rtest(deg2rad(abs(compAngDiffs)).*2);
    fprintf(fid,'\n\n\tRayleighs Test p = %0.2e, z = %0.2f\n',pval,zval);
    fprintf(fid,'\n\tCirc Mean = %0.2f',rad2deg(circ_mean(deg2rad(abs(compAngDiffs)))));
    fprintf(fid,'\n\tCirc Standard Deviation = %0.2f\n',rad2deg(circ_std(deg2rad(abs(compAngDiffs)))));
    fclose(fid);
    
end




















