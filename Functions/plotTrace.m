function plotTrace(paths)
    for p = paths'
        s = load(p{1});
        
%         v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
%         m = mkTraceMaps(s.processed.p,s.processed.trace,[]);
%         m(:,:,s.processed.splithalf.p > 0.05) = [];
        
%         gT = s.processed.trace(s.processed.splithalf.p<=0.05,:);
        gT = s.processed.trace;
        
        [isIn isMostRecent] = isInROI(s.processed.p,s.processed.roi.door);
%         
%         com = isMostRecent(1,:)|isMostRecent(3,:)|isMostRecent(5,:);
%         uni = isMostRecent(2,:)|isMostRecent(4,:)|isMostRecent(6,:);
        
        doK = [3 3];
        for part = 0:floor(length(gT(:,1))/prod(doK))

            step = 30;
            
            figure(1)
            set(gcf,'position',[50 50 1800 1200])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(gT(:,1))
                    break
                end
                subplot(doK(1),doK(2),k)
                
                v = sqrt(nansum(diff(s.processed.p,[],2).^2)).*30; 
                
                plot(s.processed.p(2,1:step:end),s.processed.p(1,1:step:end), ...
                    'color',[0.8 0.8 0.8],'linewidth',1) %[0.8 0.8 0.8]
                hold on
                
                plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))),...
                    s.processed.p(1,logical(gT(part.*prod(doK)+k,:))),...
                    'color',[0 0 0],'linestyle','none','marker','o',...
                    'markerfacecolor',[0 0 0],'markersize',2);
                
%                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&com),...
%                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&com),...
%                     'color',[0.5 0.1 1],'linestyle','none','marker','o',...
%                     'markerfacecolor',[0.5 0.1 1],'markersize',2);
% 
%                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&uni),...
%                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&uni),...
%                     'color',[0.1 1 1],'linestyle','none','marker','o',...
%                     'markerfacecolor',[0.1 1 1],'markersize',2);
                
%                 drawnow

%                 title(num2str(part.*prod(doK)+k))
                axis tight
                axis equal
                axis off
                
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/TracePlots_Undifferentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff');
            outP = ['Plots/TracePlots_Undifferentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'pdf');
            close all
            drawnow
            
            close all
            figure(1)
            set(gcf,'position',[50 50 1400 900])
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(gT(:,1))
                    break
                end
                subplot(doK(1),doK(2),k)
                
                plot(s.processed.p(2,1:step:end),s.processed.p(1,1:step:end), ...
                    'color',[0.0 0.0 0.0],'linewidth',1) %[0.8 0.8 0.8]
                hold on
                
%                 plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))),...
%                     s.processed.p(1,logical(gT(part.*prod(doK)+k,:))),...
%                     'color',[0 0 0],'linestyle','none','marker','o',...
%                     'markerfacecolor',[0 0 0],'markersize',2);
                
                plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&isMostRecent(1,:)),...
                    s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&isMostRecent(1,:)),...
                    'color',[0.5 0.1 1],'linestyle','none','marker','o',...
                    'markerfacecolor',[0.5 0.1 1],'markersize',2);

                plot(s.processed.p(2,logical(gT(part.*prod(doK)+k,:))&isMostRecent(2,:)),...
                    s.processed.p(1,logical(gT(part.*prod(doK)+k,:))&isMostRecent(2,:)),...
                    'color',[0.1 1 1],'linestyle','none','marker','o',...
                    'markerfacecolor',[0.1 1 1],'markersize',2);
                
%                 drawnow

%                 title(num2str(part.*prod(doK)+k))
                axis tight
                axis equal
                axis off
                
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/TracePlots_Differentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff');
            outP = ['Plots/TracePlots_Differentiated/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'pdf');
            close all
            drawnow
        end
    end
end