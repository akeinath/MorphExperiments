function out = chooseEmbeddingExamples(clustered,doV)

    c2c_embedding = clustered.embedding;

    K = 5;
    
    % By hand embedding
    ass = nan(length(c2c_embedding(:,1)),1);
    clustCount = 0;
    keepGoing = true;
    while keepGoing
        clustCount = clustCount+1;
        figure(1)
        set(gcf,'position',[50 50 600 600])
        scatter(c2c_embedding(isnan(ass),1),c2c_embedding(isnan(ass),2),20,'k','filled');
        axis equal
        [px py] = getpts();
        
        if isempty(px)
            keepGoing = false;
        else
            best = knnsearch(c2c_embedding,[px py],'K',K);
        
            figure(2)
            set(gcf,'position',[500 50 150.*(K+1) 150])
            for i = 1:K
                subplot(1,K+1,i)
                toPlot = squarify(doV(:,:,best(i)));
                imagesc(toPlot)
                alpha(double(~isnan(toPlot)))
                axis equal
                axis off
                caxis([-1 1])
            end
            subplot(1,K+1,K+1)
%             figure(2)
%             set(gcf,'position',[700 50 150 150])
            toPlot = squarify(nanmean(doV(:,:,best),3));
            imagesc(toPlot)
            alpha(double(~isnan(toPlot)))
            axis equal
            axis off
            caxis([-1 1])
            colormap magma
        end
    end
    close all
    drawnow
    kVal = nanmax(ass);
    
    
    clear out
    [a b] = sort(ass);
    out.inds = repmat({[]},[1 kVal]);
    for i = 1:kVal
        out.inds{i} = cInds(ass==i);
    end
    [a orderClusters] = sort(cellfun(@length,out.inds),'descend');
    
    out.inds = out.inds(orderClusters);
    newAss = nan(size(ass));
    for i = 1:kVal
        newAss(ass==orderClusters(i)) = i;
    end
    out.assignment = newAss;
    out.params = params;
    out.embedding = c2c_embedding;
    
    tmp = toc;
    fprintf('  %0.3fs.',tmp);
end