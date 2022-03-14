function helpPlotMaps(m,outP)
    doK = [8 4];

    toPlot = [];
    for i = 1:length(m(1,1,1,:))
        if i > 1
            toPlot = cat(2,toPlot,nan(size(m(:,1,:,1))));
        end
        norm = repmat(nanmax(nanmax(m(:,:,:,i),[],1),[],2),[size(m(:,:,1,i))]);
        toPlot = cat(2,toPlot,m(:,:,:,i)./norm);
%         toPlot = cat(2,toPlot,m(:,:,:,i)./1);
    end
    
    for part = 0:floor(length(toPlot(1,1,:,1))/prod(doK))

        figure(1)
        set(gcf,'position',[50 50 300.*length(m(1,1,1,:)) 700])
        for k = 1:prod(doK)
            if part.*prod(doK)+k > length(toPlot(1,1,:))
                break
            end

            subplot(doK(1),doK(2),k)
            imagesc(toPlot(:,:,part.*prod(doK)+k))
            colormap parula
            alpha(double(~isnan(toPlot(:,:,part.*prod(doK)+k))))
            axis equal
            axis off    
        end
        outP2 = [outP '_Part_' num2str(part)];
        saveFig(gcf,outP2,[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
end