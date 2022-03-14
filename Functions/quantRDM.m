function quantRDM(across,root)

    conStack = [];
    for i = 1:5
        conStack = cat(3,conStack,permute(across.RDMs.partitioned(:,:,i,i,:),[1 2 3 5 4]));
    end
    nConComps = permute(nansum(nansum(~isnan(conStack),1),2),[4 3 1 2]);
    isFullCon = nConComps==nanmax(nConComps(:));
    [x y] = ind2sub(size(nConComps),find(isFullCon));
    
    doV = nan(6,6,length(x));
    for k = 1:length(x)
        doV(:,:,k) = across.RDMs.partitioned(:,:,y(k),y(k),x(k));
    end
    
    
%     doV = across.simXcon;
%     nConComps = permute(nansum(nansum(~isnan(doV),1),2),[4 3 1 2]);
%     isFullCon = nConComps==nanmax(nConComps(:));
%     doV = doV(:,:,isFullCon);
    
%     tic
%     [contextParams contextError] = fitContextPattern(doV,1);
%     toc
%     
%     save('IndividualContextFits','contextParams','contextError');
    
%     driftParams = getSameEnvDrift(across);
    
%     tic
%     predictContext(doV,root);
%     toc
    
    load('IndividualContextFits.mat');
    conMAE = nan(size(nConComps));
    for k = 1:length(x)
        conMAE(x(k),y(k)) = contextError(k);
    end
%     
%     tic
%     [blah contextPolyError] = fitDriftPattern(doV,2);
%     toc
%     
    
%     save('IndividualContextFits_Poly','contextPolyError');
    
%     
%     save('IndividualContextFits','contextParams','contextError');

    tic
    [out linMAE] = help_predCellDOL(doV,0:4);
    toc
    
    tic
    [out contextCrossvalError] = help_predCellDOS(doV,1);
    toc
    
    sVals = doV;
    for k = 1:length(doV(1,1,:))
        o = randperm(length(doV(:,1,1)));
        tmp = squarify(doV(o,o,k));
        tmp(isnan(doV(:,:,1))) = nan;
        sVals(:,:,k) = tmp;
    end
    
    tic
    [out slinMAE] = help_predCellDOL(sVals,[0:4]);
    toc
    
    tic
    [shuffle_params shuffleContextCrossvalError] = help_predCellDOS(sVals,1);
    toc
    
    clc
    for i = 1:50
        fprintf([ '\n\t' num2str(i)])
        tic
%         [shuffle_params shuffleContextCrossvalError(:,i)] = fitContextPattern(sVals,1);

        [blah shuffleContextPolyError(:,i)] = fitDriftPattern(doV,2);
        save('IndividualContextFits_Null_Poly','shuffleContextPolyError');
        toc
    end
    
    
    
    ae = [];
    for i = 1:4
        ae = [ae linMAE(:,i) slinMAE(:,i)];
    end
%     save('NewShuffles_DimsSeparate')
%     
%     load('NewShuffles')


    load('IndividualContextFits');
    load('IndividualContextFits_Null');
    load('IndividualContextFits_Null_Poly');
    load('IndividualContextFits_Poly');

    ae = [{log10(shuffleContextPolyError(:))} {log10(contextPolyError(:))} ...
        {log10(shuffleContextCrossvalError(:))} {log10(contextError(:))}];
    
    figure
    set(gcf,'position',[50 50 175 225])
    mkBow(ae)
    ylabel('log10(Mean squared error)')
    set(gca,'xticklabel',[{'Shuffled Poly'} {'Poly'} {'Shuffled Sig'} {'Sig'}])
    saveFig(gcf,[root '/MAE_log'],[{'tiff'} {'pdf'}])    
    
    pairwiseNonparametrics(ae,[root '/MAE_']);

    ae = [ae contextCrossvalError shuffleContextCrossvalError];
    
    close all
    figure(1)
    set(gcf,'position',[50 50 300 225])
    mkBow(log10(ae))
    ylabel('log_1_0(MAE)')
    xlabel('Order of DOL fit')
    saveFig(gcf,[root '/MAE_log'],[{'tiff'} {'pdf'}])    
%     a = tsne(contextParams);
    


    [a b] = mat2lag(doV);

    
    errorThreshold = inf;
    clusteredContext = clustParams(contextParams(contextError<errorThreshold,:),[],'mahalanobis');
    plotClustProps(clusteredContext,doV(:,:,contextError<errorThreshold),[root '/Clustering_Stability/']);
%     load('Clustered_Context');
%     save('Clustered_Context','clusteredContext');
    getEmbeddingExamples(clusteredContext,doV(:,:,contextError<errorThreshold),[root '/Clustering_Stability/']);
    
    load('Clustered_Context.mat')
    plotClustProps(clusteredContext,doV,[root '/Clustering_Stability/']);

    figure(1)
    set(gcf,'position',[50 50 1500 600])
    for i = 1:5
        subplot(2,3,i)
        doColors = v2rgb(contextParams(:,i));
        scatter(clusteredContext.embedding(:,1), ...
            clusteredContext.embedding(:,2),10,doColors,'filled');
        axis equal
        axis square
    end
    subplot(2,3,6)
    doColors = v2rgb(log10(contextError));
    scatter(clusteredContext.embedding(:,1), ...
        clusteredContext.embedding(:,2),10,doColors,'filled');
    axis equal
    axis square
    saveFig(gcf,[root '/Drift_UMAP_ParameterProjection'],[{'pdf'} {'tiff'}])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    doC = v2rgb(log10(mMAE));
    scatter(driftParams(:,1),driftParams(:,2),10,doC,'filled')
    
    scatter3(driftParams(:,1),driftParams(:,2),driftParams(:,3))
    
    figure
    set(gcf,'position',[50 50 300.*length(driftParams(1,:)) 300])
    mMAE = nanmean(conMAE,2);
    rvals = nan(4,1);
    pvals = nan(4,1);
    for i = 1:length(driftParams(1,:))

        isGood = ~isnan(driftParams(:,i)) & ~isnan(mMAE);
        
        if ~any(isGood)
            continue
        end

        subplot(1,length(driftParams(1,:)),i)
        scatter(sqrt(mMAE),driftParams(:,i),10,[0.9 0.6 0.2],'filled')
%         set(gca,'xlim',[-3 1],'ylim',[-1 3])
        lsline
        xlabel('Context fit error')
        ylabel(['Drift fit parameter' num2str(i)])
        axis square
        title(['Lag ' num2str(i)])
        [rvals(i) pvals(i)] = corr(mMAE(isGood),driftParams(isGood,i),'type','spearman');

        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.05;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.15;
        text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals(i) pvals(i) nansum(isGood)]), ...
            'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','right');
    end
    saveFig(gcf,[root 'ContextError_vs_DriftParams'],[{'pdf'} {'tiff'}])
    
    
    %%%%%% Crossvalidated context versus stability
    numActiveSessions = nansum(~isnan(across.mfr),2);
    doMFR = nanmean(across.mfr,2);    
    
    
    lagVals = getConVsStab_Matched(conMAE,across.simXseq,across.mfr,root);
    
    doLag = 4;
    kruskalwallis(cat(1,lagVals{:,doLag}),[ones(length(lagVals{1,doLag}),1); ones(length(lagVals{2,doLag}),1).*2; ...
        ones(length(lagVals{3,doLag}),1).*3; ones(length(lagVals{4,doLag}),1).*4])
    
    for doLag = 1:4
        tmp = lagVals(:,doLag);
        
        rvals = nan(length(tmp));
        pvals = nan(length(tmp));
        statsval = nan(length(tmp));
        for gi = 1:length(tmp)
            for gj = gi+1:length(tmp)
                [a b c] = ranksum(tmp{gi},tmp{gj});
                pvals(gi,gj) = a;
                statsval(gi,gj) = c.zval;
                rval(gi,gj) = c.ranksum;
            end
        end
        
        % 
        pvals < 0.05./(6)
        pvals < 0.01./(6)
        pvals < 0.001./(6)
    end
    
    
    lagVals = getConVsStab(conMAE,across.simXseq);
    
    figure
    set(gcf,'position',[50 50 1200 300])
    rvals = nan(4,1);
    pvals = nan(4,1);
    for i = 1:4

        subplot(1,4,i)
        scatter(log10(lagVals{i}(:,1)),(lagVals{i}(:,2)),10,[0.9 0.6 0.2],'filled')
        set(gca,'xlim',[-5 0],'ylim',[-0.5 1]);
        hold on
        plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
        lsline
        xlabel('Context fit error log_1_0(MAE)')
        ylabel('Stability atanh(r)')
        axis square
        title(['Lag ' num2str(i)])
        [rvals(i) pvals(i)] = corr(lagVals{i}(:,1),lagVals{i}(:,2),'type','spearman');

        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.05;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.15;
        text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals(i) pvals(i) length(lagVals{i}(:,1))]), ...
            'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','right');
    end    
    saveFig(gcf,[root 'ContextError_vs_Stability'],[{'pdf'} {'tiff'}])
    
    %%%%%% Crossvalidated stability versus stability
    
    lagVals = getCrossvalStab(across.simXseq);
    arvals = nan(4,4);
    apvals = nan(4,4);
    for i = 1:4
        for j = 1:4
            if isempty(lagVals{i,j})
                continue
            end
            [arvals(i,j) apvals(i,j)] = corr(lagVals{i,j}(:,1),lagVals{i,j}(:,2),'type','spearman');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    figure
    set(gcf,'position',[50 50 1200 1500])
    toPlot = repmat({[]},[1 4]);
    for useSeq = 1:5
    
        mMAE = conMAE(:,useSeq);
        
        tmp = across.simXseq;
        tmp(useSeq,:,:) = nan;
        tmp(:,useSeq,:) = nan;
        driftXlag = mat2lag(tmp);

        rvals = nan(4,1);
        pvals = nan(4,1);
        for i = 1:4
            
            isGood = ~isnan(driftXlag(:,i)) & ~isnan(mMAE);
            if ~any(isGood)
                continue
            end
            
            toPlot{i} = cat(1,toPlot{i},[mMAE(isGood) driftXlag(isGood,i)]);
            
            subplot(5,4,(useSeq-1).*4+i)
            scatter(log10(mMAE),atanh(driftXlag(:,i)),10,[0.9 0.6 0.2],'filled')
            set(gca,'xlim',[-3 1],'ylim',[-1 3])
            lsline
            xlabel('Context fit error log_1_0(MAE)')
            ylabel('Stability atanh(r)')
            axis square
            title(['Lag ' num2str(i)])
            [rvals(i) pvals(i)] = corr(mMAE(isGood),driftXlag(isGood,i),'type','spearman');
            
            textX = get(gca,'xlim');
            textX = textX(2) - [textX(2)-textX(1)].*0.05;
            textY = get(gca,'ylim');
            textY = textY(2) - [textY(2)-textY(1)].*0.15;
            text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals(i) pvals(i) nansum(isGood)]), ...
                'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','right');
        end
    end
    saveFig(gcf,[root 'ContextMAE_vs_DriftLag_LeaveEachOut'],[{'pdf'} {'tiff'}])
    
    figure
    set(gcf,'position',[50 50 1200 300])
    mMAE = nanmean(conMAE,2);
    driftXlag = mat2lag(across.simXseq);
    rvals = nan(4,1);
    pvals = nan(4,1);
    for i = 1:4

        isGood = ~isnan(driftXlag(:,i)) & ~isnan(mMAE);
        if ~any(isGood)
            continue
        end

        subplot(1,4,i)
        scatter(log10(toPlot{i}(:,1)),atanh(toPlot{i}(:,2)),10,[0.9 0.6 0.2],'filled')
        set(gca,'xlim',[-3 1],'ylim',[-1 3])
        lsline
        xlabel('Context fit error log_1_0(MAE)')
        ylabel('Stability atanh(r)')
        axis square
        title(['Lag ' num2str(i)])
        [rvals(i) pvals(i)] = corr(toPlot{i}(:,1),toPlot{i}(:,2),'type','spearman');

        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.05;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.15;
        text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals(i) pvals(i) nansum(isGood)]), ...
            'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','right');
    end
    saveFig(gcf,[root 'ContextMAE_vs_DriftLag_LeaveEachOutAndCombine'],[{'pdf'} {'tiff'}])
end