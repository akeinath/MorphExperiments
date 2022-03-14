function batchPlotPaths(uP,root)
    tic    

    figure
    set(gcf,'position',[50 50 8.*150 ceil(length(uP)./8).*150],'color','w')
    fprintf(['\tPlotting paths... '])
    for si = 1:length(uP)
        subplot(ceil(length(uP)./8),8,si)
        plot(uP{si}(2,1:5:end),uP{si}(1,1:5:end),'color','k','linewidth',0.5)
        axis square
        set(gca,'xlim',[-2 40],'ylim',[-2 40])
        axis off
        text(19,44,['Session ' num2str(si)],'fontname','arial','fontsize',9,...
            'horizontalalignment','center')
    end  
    saveFig(gcf,root,[{'tiff'} {'pdf'}])
    close all
    drawnow
    
    tmp = toc;
    fprintf('  %0.3fs.',tmp);
end