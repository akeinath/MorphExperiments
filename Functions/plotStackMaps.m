function plotStackMaps(m)

    doM = [4 8];
    figure
    set(gcf,'position',[50 50 doM(2).*150 doM(1).*150])
    for i = 1:doM(1)
        
        doI = (i-1).*doM(2)+1:(i).*doM(2);
        
        tcm = [];
        for j = doI
            if isempty(tcm)
                tcm = m(:,:,j)./nanmax(nanmax(m(:,:,j)));
            else
                tcm = cat(2,tcm,nan(length(tcm(:,1)),1),m(:,:,j)./nanmax(nanmax(m(:,:,j))));
            end
        end
        
        subplot(doM(1),1,i)
        imagesc(tcm)
        colormap parula
        alpha(double(~isnan(tcm)))
        for j = 1:doM(2)
            text((j-1).*(length(m(1,:,1))+1)+2.5,2.5,num2str((i-1).*doM(2)+j),...
            'fontname','arial','fontsize',9,'horizontalalignment','center');
        end
        axis equal
        axis off
    end
end