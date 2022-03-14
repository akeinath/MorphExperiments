function rdmClusterExplosion(across,root)
    
    [blah vals] = mat2lag(across.simXseq);
    doExclude = any(isnan(vals),2);
    cInds1 = find(~doExclude);
    
    [blah vals] = mat2lag(across.simXcon);
    doExclude = any(isnan(vals),2);
    cInds2 = find(~doExclude);

    doCells = intersect(cInds1,cInds2);
    
    context_clustering = seqSimClustering(across.simXcon(:,:,doCells), ...
        'raw',6,[root '/Clustering_Context/']);

    stability_clustering = seqSimClustering(across.simXseq(:,:,doCells), ...
        'lag',6,[root '/Clustering_Stability/']);
    
    
    predictContext(across.simXcon(:,:,doCells),[root '/Predicting_Context_RDM_DOS']);
    predictContext(across.simXseq(:,:,doCells),[root '/Predicting_Drift_RDM_DOS']);
    
    comparePredictions(across,doCells,[root '/ComparePredictions'])
    
    transCounts = nan(length(stability_clustering.inds),length(context_clustering.inds));
    for i = 1:length(stability_clustering.inds)
        for j = 1:length(context_clustering.inds)
            transCounts(i,j) = length(intersect(stability_clustering.inds{i},...
                context_clustering.inds{j}));
        end
    end
    
    figure()
    set(gcf,'position',[50 50 300 600])
    alluvialflow(transCounts,num2cell(1:length(stability_clustering.inds)), ...
        num2cell(1:length(context_clustering.inds)))
    saveFig(gcf,[root '/Alluvial'],[{'tiff'} {'pdf'} {'jpeg'}])
    
%     save('MatlabData/BatchClustering_Context.mat','-struct','context_clustering','-v7.3');
    save('MatlabData/BatchClustering_Stability.mat','-struct','stability_clustering','-v7.3');
    
    context_clustering = load('MatlabData/BatchClustering_Context.mat');
    stability_clustering = load('MatlabData/BatchClustering_Stability.mat');
%     
%     [rval pval] = corr(stability_clustering.params(:,1),context_clustering.error);

    figure
    set(gcf,'position',[50 50 600 300])
    subplot(1,2,1)
    scatter(stability_clustering.rotEmbedding(:,1),context_clustering.error,15,'k','filled');
    ylabel('DOS fit mean absolute error')
    xlabel('Drift fit PC1')
    lsline
    axis square
    subplot(1,2,2)
    scatter(stability_clustering.rotEmbedding(:,2),context_clustering.error,15,'k','filled');
    ylabel('DOS fit mean absolute error')
    xlabel('Drift fit PC2')
    axis square
    lsline
    
    
    [rval pval] = corr(stability_clustering.rotEmbedding(:,1),context_clustering.error);
    [rval pval] = corr(stability_clustering.rotEmbedding(:,2),context_clustering.error);
end
    