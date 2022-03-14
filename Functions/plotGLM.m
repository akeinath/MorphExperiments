function pVal = plotGLM_2D(across,root)
    act = sqrt(across.class.d3.tr2.*across.class.d3.r2penalty);
    shuf = sqrt(across.class.d3.shuffle_tr2.*across.class.d3.shuffle_r2penalty);
    
    isGood = ~any(isnan(act),2) & ~any(isnan(shuf),2);
    act(~isGood,:) = [];
    shuf(~isGood,:) = [];
    
    actD = mahal(act,shuf);
    shufD = mahal(shuf,shuf);
    pVal = 1-(nanmean(bsxfun(@gt,actD,shufD'),2));
    shufPVal = 1-(nanmean(bsxfun(@gt,shufD,shufD'),2));
    
    cm = flipud(hot(101)).*0.75;
    toMap = 0.05-pVal;
    toMap(toMap<0) = 0;
    toMap = round(toMap.*2000);
    
    toMapShuf = 0.05-shufPVal;
    toMapShuf(toMapShuf<0) = 0;
    toMapShuf = round(toMapShuf.*2000);
    
    figure
    set(gcf,'position',[50 50 600 300])
    subplot(1,2,1)
    scatter(act(:,1),act(:,2),3,cm(toMap+1,:))
%     scatter(act(:,1),act(:,2),3,'markeredgecolor',[0.6 0.6 0.9])
    set(gca,'xlim',[0 1],'ylim',[0 1])
    xlabel('Drift-attributable r','fontname','arial')
    ylabel('Context-attributable r','fontname','arial')
    axis square
    subplot(1,2,2)
    scatter(shuf(:,1),shuf(:,2),3,cm(toMapShuf+1,:))
    set(gca,'xlim',[0 1],'ylim',[0 1])
    xlabel('Drift-attributable r','fontname','arial')
    ylabel('Context-attributable r','fontname','arial')
    axis square
    saveFig(gcf,[root '/Context_vs_Drift_r2'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 50 250 250])
    cumHist(pVal,[0:0.01:1]);
    hold on
    plot([0 1],[0 1],'linestyle','--','color','k');
    axis square
    saveFig(gcf,[root '/Mahal_pVal'],[{'pdf'} {'tiff'}]);
end