function plotColoredEmbedding(c2c_embedding,doColors,root)   
    if nargin < 2 || isempty(doColors)
        doColors = ones(length(c2c_embedding(:,1)),3).*0.5;
    end

    figure(1)
    set(gcf,'position',[50 50 400 400]);
    scatter(c2c_embedding(:,1)-1,c2c_embedding(:,2),30,doColors,'filled')
%     set(gca,'xlim',[-2 0],'ylim',[-1 1]);
    hold on
    if length(c2c_embedding(1,:))>2
        set(gcf,'position',[50 50 1000 400]);
        scatter(c2c_embedding(:,1)+1,c2c_embedding(:,3),30,doColors,'filled')
        scatter(c2c_embedding(:,2)+3,c2c_embedding(:,3),30,doColors,'filled')
        set(gca,'xlim',[-2 4],'ylim',[-1 1]);
        axis equal
        plot([0 0],get(gca,'ylim'),'linestyle','--','color','k')
        plot([2 2],get(gca,'ylim'),'linestyle','--','color','k')
    end
%     h(1) = plot([0]',[-inf]','markersize',10,'color',hulkcm([1 0],[0 1]), ...
%         'markerfacecolor',hulkcm([1 0],[0 1]),'linestyle','none','marker','o');
%     h(2) = plot([0]',[-inf]','markersize',10,'color',hulkcm([0 1],[0 1]), ...
%         'markerfacecolor',hulkcm([0 1],[0 1]),'linestyle','none','marker','o');
%     legend(h,[{'Drift-loading'} {'Context-loading'}],'location','south','fontname','arial',...
%         'fontsize',9);
    saveFig(gcf,[root '/Cell2Cell_Isomap/Time_vs_Context'],[{'tiff'} {'pdf'} {'jpeg'}]);
end