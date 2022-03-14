function out = seqSimClustering(seqSim,lagVsraw,kVal,root)

    close all
    drawnow

    if isempty(kVal)
        kVal = 2;
    end

    if isempty(lagVsraw)
        lagVsraw = {'raw'};
    end
    
    [vals] = mat2lag(seqSim);
    cInds = find(~any(isnan(vals),2));
    
    out = clustParams(seqSim(:,:,cInds),cInds,lagVsraw,kVal,'minkowski');
    
    figure(2)
    [a b] = sort(out.assignment);
    set(gcf,'position',[50 50 500 500]);
    subplot(1,5,1:3)
    imagesc(vals(b,:))
    linLocs = cumsum(cellfun(@length,out.inds));
    hold on
    plot(repmat(get(gca,'xlim'),[length(linLocs)-1 1])', ...
        0.5+[linLocs(1:end-1)' linLocs(1:end-1)']','color','k',...
        'linewidth',0.5)
    caxis([-1 1])
    colorbar()
    subplot(1,5,4:5)
    plot(vals(b,:),flipud([1:length(vals(:,1))]'),'color','k')
    plot(nanmean(vals(b,:),2),flipud([1:length(vals(:,1))]'),'color','k')
    hold on
    plot([0 0],[1 length(vals(:,1))],'linestyle','--','color',[0.5 0.5 0.5]);
    set(gca,'xlim',[-1 1],'ylim',[0 length(vals(:,1))+1])
    drawnow
    saveFig(gcf,[root '/List'],[{'pdf'} {'tiff'} {'jpeg'}]);

    if isfield(out,'embedding')
        c2c_embedding = out.embedding;
    else
%         opt.dims = 2;
%         opt.display = false;
%         embedding = IsoMap(squareform(pdist(out.params,'mahalanobis')),'k',30,opt);
        c2c_embedding = out.params; %embedding.coords{1}';
    end
    
    clear h
    figure(1)
    hold off
    set(gcf,'position',[50 50 500 500])
    colors =(jet(length(out.inds)+2));
%     colors(6,:) = [];

    toPlot = [c2c_embedding*pca(c2c_embedding)];
    out.rotEmbedding = toPlot;
%     doColors = v2rgb(toPlot(:,2));
%     scatter(c2c_embedding(:,1),c2c_embedding(:,2),30,doColors,'filled');
        
    for i = 1:length(out.inds)
        h(i) = scatter(c2c_embedding(out.assignment==i,1), ...
            c2c_embedding(out.assignment==i,2),30,colors(i+1,:),'filled');
        hold on
    end
    for i = 1:length(out.inds)
        text(nanmedian(c2c_embedding(out.assignment==i,1)), ...
            nanmedian(c2c_embedding(out.assignment==i,2)),num2str(i),...
            'fontname','arial','color','k','fontweight','bold','fontsize',12, ...
            'horizontalalignment','center','verticalalignment','middle');
    end
    saveFig(gcf,[root '/IsomapEmbedding'],[{'pdf'} {'tiff'} {'jpeg'}]);

    figure(3)
    set(gcf,'position',[50 50 400 200.*length(out.inds)])
    for i = 1:length(out.inds)
        subplot(length(out.inds),2,((i-1).*2)+1)
        imagesc(squarify(nanmean(seqSim(:,:,out.inds{i}),3)));
        caxis([0.15 0.7])
        alpha(double(~eye(size(seqSim(:,:,i)))));
        axis equal
        axis off
        subplot(length(out.inds),2,((i-1).*2)+2)
        mkLine(mat2lag(seqSim(:,:,out.inds{i})))
        hold on
        set(gca,'ylim',[-0.2 0.9],'xlim',[0.75 length(seqSim(:,1,1))-1+0.25])
        plot(get(gca,'xlim'),[0 0],'linestyle','--','color',[0.5 0.5 0.5]);
    end
    saveFig(gcf,[root '/CollapsedRDMView'],[{'pdf'} {'tiff'} {'jpeg'}]);
end








