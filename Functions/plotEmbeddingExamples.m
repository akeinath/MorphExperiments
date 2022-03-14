function plotEmbeddingExamples(clustered,ex)

    c2c_embedding = clustered.embedding;

   
%     figure(1)
%     set(gcf,'position',[50 50 600 600])
    scatter(c2c_embedding(:,1),c2c_embedding(:,2),10,[0.75 0.75 0.75],'filled');
    axis equal
    hold on
    scatter(c2c_embedding(ex,1),c2c_embedding(ex,2),10,'k','filled');
    for k = 1:length(ex)
        text(c2c_embedding(ex(k),1),c2c_embedding(ex(k),2)+0.25,num2str(k),'fontname','arial',...
            'fontsize',5,'horizontalalignment','center')
    end

end