function plotWeightClass(r2penalty,tr2,outP,firstIsOverall)

    if nargin<4 || isempty(firstIsOverall)
        firstIsOverall = false;
    end

    colors = bsxfun(@times,[1 0 0],r2penalty(:,1));
    colors = colors+bsxfun(@times,[0 0 1],r2penalty(:,2));
    colors = colors+bsxfun(@times,[0 1 0],r2penalty(:,3));
    colors(colors<0) = 0;
    colors = colors.^(1./4);
    colors = [colors./nanmax(colors(:))];

    if firstIsOverall
        overallR2Penalty = r2penalty(1,:);
        overallTR2 = tr2(1,:);
        overallColor = colors(1,:);
        r2penalty(1,:) = [];
        tr2(1,:) = [];
        colors(1,:) = [];
    end

    [blah best] = nanmax(tr2(~isnan(tr2)));
    
    %%%%%%%%% Plot 3D Gif
    
% % %     figure
% % %     set(gcf,'position',[50 50 600 600])
% % %     hold on
% % %     scatter3(r2penalty(:,1),r2penalty(:,2),r2penalty(:,3),(200.*tr2(~isnan(tr2))),colors,'filled');
% % %     
% % % %     tmpTR = tr2(~isnan(tr2));
% % % %     scatter3(r2penalty(best,1),r2penalty(best,2),r2penalty(best,3),(200.*tmpTR(best))./2,'k',...
% % % %         'filled');
% % %     
% % %     if firstIsOverall
% % %         scatter3(overallR2Penalty(:,1),overallR2Penalty(:,2),overallR2Penalty(:,3), ...
% % %             (500.*overallTR2),overallColor,'filled','marker','p');
% % %         scatter3(overallR2Penalty(:,1),overallR2Penalty(:,2),overallR2Penalty(:,3), ...
% % %             (500.*overallTR2),'k','linewidth',1.5,'marker','p');
% % %         text(0.8,0.8,0.72,['Pop r^2 = ' num2str(overallTR2)],...
% % %             'horizontalalignment','center','color',[0.25 0.25 0.25],'fontsize',12,'fontweight','bold')  
% % %     end
% % %     
% % % %     text(r2penalty(best,1),r2penalty(best,2),r2penalty(best,3)+0.06,['r^2 = ' num2str(nanmax(tr2))],...
% % % %         'horizontalalignment','center','color',[0.25 0.25 0.25])
% % %     text(0.8,0.8,0.8,['Max Cell r^2 = ' num2str(nanmax(tr2))],...
% % %         'horizontalalignment','center','color',[0.25 0.25 0.25],'fontsize',12,'fontweight','bold')  
% % %     
% % % 
% % %     zlabel('Geometry RDM dropout \Deltar^2 (%)');
% % %     ylabel('Attractor RDM dropout \Deltar^2 (%)');
% % %     xlabel('Time RDM dropout \Deltar^2 (%)');
% % %     axis equal
% % %     set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',[0 1])
% % %     grid on
% % %     clear frames
% % %     for i = 1:1:360
% % %         view(i+0.5,30)
% % %         drawnow
% % %         frames(i) = getframe(gcf);
% % % 
% % %         im = frame2im(frames(i)); 
% % %         [imind,cm] = rgb2ind(im,256); 
% % %         if i == 1
% % %             checkP(outP);
% % %             imwrite(imind,cm,outP,'gif', 'Loopcount',inf,'DelayTime',0.025); 
% % %         else
% % %             imwrite(imind,cm,outP,'gif','DelayTime',0.025,'WriteMode','append'); 
% % %         end
% % %     end
    
    %%%%%%%%% Plot Ternary Plot

    figure
    set(gcf,'position',[50 50 400 400],'color','w') 
    plotTernary(r2penalty,(100.*tr2(~isnan(tr2))),colors,[{'Time'} {'Attractor'} {'Geometry'}]);
    hold on
    text(-1,0.9,['Max Cell r^2 = ' num2str(round(1000.*nanmax(tr2))./1000)],...
        'horizontalalignment','left','color','k','fontsize',9)
    text(-1,0.8,['n = ' num2str(nansum(~isnan(tr2)))],...
        'horizontalalignment','left','color','k','fontsize',9)
    if firstIsOverall
        h = plotTernary(overallR2Penalty,(300.*overallTR2),overallColor,[],false);
        set(h,'marker','p','markerfacecolor','w')
        text(-1,1,['Population r^2 = ' num2str(round(1000.*nanmax(overallTR2))./1000)],...
            'horizontalalignment','left','color','k','fontsize',9) 
    end
    tOutP = [outP(1:end-4) '_TernaryPlot'];
    saveFig(gcf,tOutP,[{'pdf'} {'tiff'}]);
    
%     figure
%     set(gcf,'position',[50 50 400 400],'color','k') 
%     plotTernaryHeatmap(r2penalty,[{'Time'} {'Attractor'} {'Geometry'}]);
%     tOutP = [outP(1:end-4) '_TernaryHeatmap'];
%     saveFig(gcf,tOutP,[{'pdf'} {'tiff'}]);
end




























