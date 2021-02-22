function valVsRDMFits(meanSplitHalfCorr,fits,minSamples)
        %%%% Correlation with firing properties
        %%% SHC vs Attractor
        colors = bsxfun(@times,[0 0.25 1],-fits(fits(:,end)>minSamples,2));
        colors(colors<0) = 0;
        colors = colors.^(1./4);
        colors = [colors./nanmax(colors(:))];
        figure
        set(gcf,'position',[50 50 300 600])
        subplot(2,1,1)
        scatter(meanSplitHalfCorr(fits(:,end)>minSamples),-fits(fits(:,end)>minSamples,2), ...
            30,colors,'filled')
        hold on
        plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
        h = lsline;
        set(h,'color',[0 0.75 0.75],'linewidth',2)
        ylabel('Attractor RDM Fit (Kendall''s \tau)')
        xlabel('Value')
        [rval pval] = corr(meanSplitHalfCorr(fits(:,end)>minSamples), ...
            -fits(fits(:,end)>minSamples,2));
        set(gca,'ylim',[-0.2 0.8])
        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.9;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.1;
        text(textX,textY,sprintf(['r = %0.3f\np = %0.4f'],[rval pval]), ...
            'fontweight','normal','fontname','arial','fontsize',10);
        axis square
        
        %%% SHC vs Time
        colors = bsxfun(@times,[1 0 0],-fits(fits(:,end)>minSamples,1));
        colors(colors<0) = 0;
        colors = colors.^(1./4);
        colors = [colors./nanmax(colors(:))];
        subplot(2,1,2)
        scatter(meanSplitHalfCorr(fits(:,end)>minSamples),-fits(fits(:,end)>minSamples,1), ...
            30,colors,'filled')
        hold on
        plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
        h = lsline;
        set(h,'color',[0.75 0 0],'linewidth',2)
        ylabel('Time RDM Fit (Kendall''s \tau)')
        xlabel('Value')
        [rval pval] = corr(meanSplitHalfCorr(fits(:,end)>minSamples), ...
            -fits(fits(:,end)>minSamples,1));
        set(gca,'ylim',[-0.2 0.8])
        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.9;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.1;
        text(textX,textY,sprintf(['r = %0.3f\np = %0.4f'],[rval pval]), ...
            'fontweight','normal','fontname','arial','fontsize',10);
        axis square
end