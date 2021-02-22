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
    compAngDiffs = [];
    sigParams = [];
    
    shuffle_r2penalty = [];
    shuffle_tr2 = [];
    
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
        uclass = nan([length(alignMap{1}(:,1)) length(sessions)]);
        uIsPC = false([length(alignMap{1}(:,1)) length(sessions)]);
        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        
        tmp = size(SFPs{1}(1,:,:));
        if doAnatomy
            usfps = nan([tmp(2:3) length(alignMap{1}(:,1)) length(sessions)]);
        end
        for si = 1:length(sessions)
            
            umfr(alignMap{1}(:,si)~=0,si) = MFRs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            upfs(alignMap{1}(:,si)~=0,si) = PFSs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            usic(alignMap{1}(:,si)~=0,si) = SICs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            ushc(alignMap{1}(:,si)~=0,si) = SHCs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            uIsPC(alignMap{1}(:,si)~=0,si) = isPC{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
            
            if doAnatomy
                SFPs{si} = permute(SFPs{si},[2 3 1]);
                usfps(:,:,alignMap{1}(:,si)~=0,si) = SFPs{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
            end
        end
        
        if doAnatomy
            msfps = nanmean(usfps,4);
        end
      
        
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
            
                
                xc = pvxcorr3(tmp1,tmp2,[1 1 1],30);
                a = 2;
                b = 2;
                c = 2;
                ivals = xcorr3transform(tmp1,tmp2,[a-ceil(length(xc(:,1,1))./2) ...
                    b-ceil(length(xc(1,:,1))./2) c-ceil(length(xc(1,1,:))./2)]);
                sim(si,sj,:) = ivals;

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
                [root '/' mouseID '/SequenceAnalyses/Sequence_' num2str(i)]);
            
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
      
        [r2penalty tr2] = classifyCellRDMswInteract(sim,-RDMs,root);


        valVsClass(nanmedian(ushc,2),r2penalty,tr2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}])
        saveFig(gcf,[root '/Classification/SHC_vs_tr2'],[{'pdf'} {'tiff'}])
        
        allR2Penalty = [allR2Penalty; r2penalty];
        allTR2 = [allTR2; tr2];
        
        
        [ta tb] = classifyCellRDMswInteract(sim,-RDMs,root,true);
        
        shuffle_r2penalty = [shuffle_r2penalty; ta];
        shuffle_tr2 = [shuffle_tr2; tb];
        
        allCellProps = [allCellProps; nanmedian(ushc,2) nanmedian(usic,2) nanmedian(umfr,2)];

        [a b] = nanmax(r2penalty,[],2);
        b(all(isnan(r2penalty),2)) = 0;
        
        uclass = bsxfun(@times,~isnan(umfr),b);
        uclass(isnan(umfr)) = nan;
        
        tmp = nan(4,length(uclass(1,:)));
        
        for j = 1:length(uclass(1,:))
            for i = 1:4
                tmp(i,j) = nansum(uclass(:,j)==(i-1))./nansum(~isnan(uclass(:,j)));
            end
        end
        
        allClassXTime = [allClassXTime {tmp}];
        
        
    end
    
    slashInds = find(ismember(paths{1},'/'));
    root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
    
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
    
    %%%%%%%%%%% Ranked factor plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ti = [0 0.3]
    
        
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
            allStats = [allStats [pval stats.ranksum]'];
        end
        
        figure
        set(gcf,'position',[50 50 225 300],'color','w')
        mkBow(toPlot,[{'1st'} {'2nd'} {'3rd'}]);
        xlabel('Factor rank')
        ylabel('Attributable r^2')
        for i = 1:3
            text(i,1.12,sprintf('Z=%.1e\np=%.1e',abs(allStats(2,i)),allStats(1,i)),...
                'horizontalalignment','center','color','k','fontsize',7)
        end
        set(gca,'outerposition',[0 0 1 0.9])
        saveFig(gcf,[root '/Overall/Ranked_factor_tr2_threshold_' num2str(ti)],[{'pdf'} {'tiff'}])
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%
    
    for ti = [0 0.3]
        plotWeightClass(allR2Penalty(~isnan(allTR2)&allTR2>ti,:),allTR2(allTR2>ti),...
            [root '/Overall/Dropout_withInteractions_r2_threshold_' num2str(ti) '.gif'],false);
    end

    valVsClass(allCellProps(:,1),allR2Penalty,allTR2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}])
    saveFig(gcf,[root '/Overall/shc_vs_tr2'],[{'pdf'} {'tiff'}])
    
    valVsClass(allCellProps(:,3),allR2Penalty,allTR2, ...
            [{'Time Cells'} {'Attractor Cells'} {'Geometry Cells'}])
    saveFig(gcf,[root '/Overall/mfr_vs_tr2'],[{'pdf'} {'tiff'}])

    fprintf('\n\n\t**************************************\n\n')
    fprintf('\tEmbedded Angle absolute differences:')
    fprintf('\n\t\t%0.2f',abs(compAngDiffs));
    [pval zval] = circ_rtest(deg2rad(abs(compAngDiffs)).*2);
    fprintf('\n\n\tRayleighs Test p = %0.2e, z = %0.2f\n',pval,zval)
    fprintf('\n\tCirc Mean = %0.2f',rad2deg(circ_mean(deg2rad(abs(compAngDiffs)))))
    fprintf('\n\tCirc Standard Deviation = %0.2f\n',rad2deg(circ_std(deg2rad(abs(compAngDiffs)))))
end




















