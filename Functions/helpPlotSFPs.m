function helpPlotSFPs(SFPs,outP)
    
    doK = [8 4];

%     keyColor = [0.9 0.1 0.1; 0.9 0.4 0.4; 0.9 0.7 0.7; 0.7 0.7 0.9; 0.4 0.4 0.9; 0.1 0.1 0.9].*0.5+0.5;
    keyColor 
    keyColor = [keyColor(1,:); keyColor; keyColor(end,:)];
    
    for part = 0:floor(length(SFPs(1,1,1,:))/prod(doK))
        
        figure
        set(gcf,'position',[50 50 125.*doK(2) 125.*doK(1)])
        
        for k = 1:prod(doK)
            if part.*prod(doK)+k > length(SFPs(1,1,1,:))
                break
            end
            
            subplot(doK(1),doK(2),k)
            hold on
            for j = 1:length(SFPs(1,1,:,1))
                if ~all(all(isnan(SFPs(:,:,j,part.*prod(doK)+k))))
                    tmp = regionprops(SFPs(:,:,j,part.*prod(doK)+k)>0,'ConvexHull');
                    if ~isempty(tmp)
                        plot(tmp(1).ConvexHull(:,1),tmp(1).ConvexHull(:,2),'color',keyColor(j,:));
                    end
                end
            end
            set(gca,'xlim',[0 length(SFPs(:,1,1,1))+1],...
                'ylim',[0 length(SFPs(:,1,1,1))+1])
            axis square
        end
        outP2 = [outP '_Part_' num2str(part)];
        saveFig(gcf,outP2,[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
end