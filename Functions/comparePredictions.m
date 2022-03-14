function comparePredictions(across,doCells,root)
    
    driftNull = load('Predicting_Drift_RDM_DOS');
    contextNull = load('Predicting_Context_RDM_DOS');
    driftNull = driftNull.null_MAE;
    contextNull =  contextNull.null_MAE;
    
    [out driftActual] = help_predCellDOS(across.simXseq(:,:,doCells),1);
    [out contextActual] = help_predCellDOS(across.simXcon(:,:,doCells),1);
    
    driftP = nanmean(bsxfun(@gt,driftActual,driftNull),2);
    contextP = nansum(bsxfun(@gt,contextActual,contextNull),2) ./ nansum(~isnan(contextNull),2);
    
    figure
    set(gcf,'position',[50 50 250 250])
    scatter(log10(driftActual),log10(contextActual),15,[0.9 0.6 0.2],'filled');
    axis square
    set(gca,'xlim',[-2 0],'ylim',[-2 0])
    lsline
    ylabel('Context log_1_0(MAE)')
    xlabel('Drift log_1_0(MAE)')
    saveFig(gcf,[root '/LogMAECorrelation'],[{'pdf'} {'tiff'}])
    
    
    figure
    set(gcf,'position',[50 50 250 250])
    scatter(contextParams(:,4),log10(contextActual),15,[0.9 0.6 0.2],'filled');
    axis square
    axis equal
%     set(gca,'xlim',[-2 0],'ylim',[-2 0])
    lsline
    ylabel('Context log_1_0(MAE)')
    xlabel('Drift log_1_0(MAE)')
    saveFig(gcf,[root '/DriftMainEffect_vs_ContextLogMAECorrelation'],[{'pdf'} {'tiff'}])
    
    figure
    set(gcf,'position',[50 50 250 250])
    scatter(driftParams(:,4),log10(contextActual),15,[0.9 0.6 0.2],'filled');
    axis square
    axis equal
%     set(gca,'xlim',[-2 0],'ylim',[-2 0])
    lsline
    ylabel('Context log_1_0(MAE)')
    xlabel('Drift log_1_0(MAE)')
    saveFig(gcf,[root '/DriftMainEffect_vs_ContextLogMAECorrelation'],[{'pdf'} {'tiff'}])
    
    
    [rval pval] = corr(driftParams(:,4),log10(contextActual(:,end)));
    
    [rval pval] = corr(log10(driftActual(:,1)),log10(contextActual(:,end)));
    
    [driftParams error] = fitContextPattern(across.simXseq(:,:,doCells),1);
    [contextParams error] = fitContextPattern(across.simXcon(:,:,doCells),1);
    
    
    driftClustering = clustParams(driftParams,doCells);
    plotClustProps(driftClustering,across.simXseq,[root '/Clustering_Stability/']);
    figure(1)
    set(gcf,'position',[50 50 750 750])
    for i = 1:4
        subplot(2,2,i)
        doColors = v2rgb(driftParams(:,i));
        scatter(driftClustering.embedding(:,1), ...
            driftClustering.embedding(:,2),20,doColors,'filled');
        axis equal
    end
    saveFig(gcf,[root '/Drift_UMAP_ParameterProjection'],[{'pdf'} {'tiff'}])
    
    contextClustering = clustParams(contextParams,doCells);
    plotClustProps(contextClustering,across.simXcon,[root '/Clustering_Stability/']);
    figure(1)
    set(gcf,'position',[50 50 750 750])
    for i = 1:4
        subplot(2,2,i)
        doColors = v2rgb(contextParams(:,i));
        scatter(contextClustering.embedding(:,1), ...
            contextClustering.embedding(:,2),20,doColors,'filled');
        axis equal
    end
    saveFig(gcf,[root '/Context_UMAP_ParameterProjection'],[{'pdf'} {'tiff'}])
    
    opt.dims = 2;
    opt.display = false;
    embedding = IsoMap(squareform(pdist(contextParams,'minkowski')),'k',30,opt);
    c2c_embedding = embedding.coords{1}';
    scatter(c2c_embedding(:,1),c2c_embedding(:,2))
    
    opt.dims = 2;
    opt.display = false;
    embedding = IsoMap(squareform(pdist(driftParams,'minkowski')),'k',30,opt);
    c2c_embedding = embedding.coords{1}';
    scatter(c2c_embedding(:,1),c2c_embedding(:,2))
    
    embedding = IsoMap(squareform(pdist(out.params,'mahalanobis')),'k',30,opt);
    
    scatter(driftParams(:,end),log10(contextActual(:,end)))
    
    scatter(driftParams(:,1),log10(contextActual(:,end)))
    
    scatter(driftParams(:,2),log10(contextActual(:,end)))
    

    close all
    [c2c_embedding b c d] = run_umap(driftParams,'n_neighbors',50, ...
        'n_epochs',500, ...
        'distance','euclidean',...
        'min_dist',0.2);
    
    figure(1)
    set(gcf,'position',[50 50 600 600])
    set(gca,'xlim',[nanmin(c2c_embedding(:))-5 nanmax(c2c_embedding(:))+5], ...
        'ylim',[nanmin(c2c_embedding(:))-5 nanmax(c2c_embedding(:))+5])
    scatter(c2c_embedding(:,1),c2c_embedding(:,2),20,'k','filled');
end