function pVal = plotGLM_3D(across,root)
    act = (across.class.d3.r2penalty); %.*across.class.d3.tr2
    shuf = (across.class.d3.shuffle_r2penalty); %.*across.class.d3.shuffle_tr2
    
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
    set(gcf,'position',[50 50 1000 1000])
    subplot(2,2,1)
    plotTernary(act,10,cm(toMap+1,:),[{'Time'} {'Context Group'} {'Shape'}],true)
    subplot(2,2,2)
    plotTernary(shuf,10,cm(toMap+1,:),[{'Time'} {'Context Group'} {'Shape'}],true)
    subplot(2,2,3)
    plotTernary(act.^2,10,cm(toMap+1,:),[{'Time'} {'Context Group'} {'Shape'}],true)
    subplot(2,2,4)
    plotTernary(shuf.^2,10,cm(toMap+1,:),[{'Time'} {'Context Group'} {'Shape'}],true)
    saveFig(gcf,[root '/Context_vs_Drift_r2'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 50 250 250])
    cumHist(pVal,[0:0.01:1]);
    hold on
    plot([0 1],[0 1],'linestyle','--','color','k');
    axis square
    saveFig(gcf,[root '/Mahal_pVal'],[{'pdf'} {'tiff'}]);
end