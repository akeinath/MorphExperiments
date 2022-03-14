function plotClustProps(out,sim,root)

    c2c_embedding = out.embedding;

    figure(1)
    hold off
    set(gcf,'position',[50 50 250 250])
    colors = 1-(inferno(length(out.inds)+2));
%     colors(6,:) = [];
%     toPlot = [c2c_embedding*pca(c2c_embedding)];
%     out.rotEmbedding = toPlot;
%     doColors = v2rgb(toPlot(:,2));
%     scatter(c2c_embedding(:,1),c2c_embedding(:,2),30,doColors,'filled');
        
    for i = 1:length(out.inds)
        h(i) = scatter(c2c_embedding(out.assignment==i,1), ...
            c2c_embedding(out.assignment==i,2),15,colors(i+1,:),'filled');
        hold on
    end
    for i = 1:length(out.inds)
        text(nanmedian(c2c_embedding(out.assignment==i,1)), ...
            nanmedian(c2c_embedding(out.assignment==i,2)),num2str(i),...
            'fontname','arial','color','k','fontweight','bold','fontsize',12, ...
            'horizontalalignment','center','verticalalignment','middle');
    end
    lim = nanmax(abs([get(gca,'xlim') get(gca,'ylim')]));
    set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
%     axis equal
    axis square
    saveFig(gcf,[root '/Embedding'],[{'pdf'} {'tiff'} {'jpeg'}]);

    figure(3)
    set(gcf,'position',[50 50 150.*ceil(sqrt(length(out.inds))) 150.*ceil(sqrt(length(out.inds)))])
    for i = 1:length(out.inds)
        subplot(ceil(sqrt(length(out.inds))),ceil(sqrt(length(out.inds))),i)
        imagesc(squarify(nanmean(sim(:,:,out.inds{i}),3)));
%         caxis([-0.2 0.8])
        caxis([-1 1])
        alpha(double(~eye(size(sim(:,:,i)))));
        axis equal
        axis off
        title(num2str(length(out.inds{i})))
%         subplot(length(out.inds),2,((i-1).*2)+2)
%         mkLine(mat2lag(sim(:,:,out.inds{i})))
%         hold on
%         set(gca,'ylim',[-0.2 0.9],'xlim',[0.75 length(sim(:,1,1))-1+0.25])
%         plot(get(gca,'xlim'),[0 0],'linestyle','--','color',[0.5 0.5 0.5]);
    end
    saveFig(gcf,[root '/CollapsedRDMView'],[{'pdf'} {'tiff'} {'jpeg'}]);
end