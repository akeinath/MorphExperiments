function plotgRDMs(gRDMs,root)
    
    maxDepth = nanmax(cellfun(@size,gRDMs,repmat({3},size(gRDMs))));
    figure
    set(gcf,'position',[50 50 300.*length(gRDMs) 300.*maxDepth]);
    for i = 1:length(gRDMs)
        for j = 1:length(gRDMs{i}(1,1,:))
            subplot(maxDepth,length(gRDMs),i+(j-1).*length(gRDMs))
            imagesc(gRDMs{i}(:,:,j))
            axis square
            axis off
            colormap inferno
        end
    end
    saveFig(gcf,root,[{'tiff'} {'pdf'}]);
end