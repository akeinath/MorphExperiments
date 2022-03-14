function seqSimWMatching(across, root)


    numActiveSessions = nansum(~isnan(across.mfr),2);
    driftParams = getSameEnvDrift(across);

    [a b] = mat2lag(across.simXseq);
    for i = 1:4
        [rval pval] = corr(numActiveSessions(~isnan(a(:,i))),a(~isnan(a(:,i)),i),'type','spearman');
        subplot(1,4,i)
        scatter(numActiveSessions(~isnan(a(:,i))),a(~isnan(a(:,i)),i));
    end
    
    
    doMFR = nanmean(across.mfr,2);    
    scatter(numActiveSessions,doMFR);
    
    [rval pval] = corr(numActiveSessions,doMFR,'type','spearman');
    
    [a b] = mat2lag(across.simXseq);
    for i = 1:2
        [rval pval] = corr(numActiveSessions(~isnan(driftParams(:,i))),driftParams(~isnan(driftParams(:,i)),i),'type','spearman');
        subplot(1,2,i)
        scatter(numActiveSessions(~isnan(driftParams(:,i))),driftParams(~isnan(driftParams(:,i)),i));
    end
    
%     av = mat2lag(across.simXseq);
    for i = 1:4
        train = i:i+1;
        test = [1:i-1 i+2:5];
        stl = mat2lag(across.simXseq(train,train,:));
        [blah ltl] = mat2lag(across.simXseq(test,test,:));
        
        isGood = ~isnan(stl(:,1)) & all(~isnan(ltl),2);
        stl = stl(isGood);
        ltl = ltl(isGood,:);
        
        rankVal = normRank(stl)./length(stl);
        
        toPlot = [{stl(rankVal<=0.5)}; {stl(rankVal>0.5)}];
        for j = 1:length(ltl(1,:))
            toPlot = [toPlot [{ltl(rankVal<=0.5,j)}; {ltl(rankVal>0.5,j)}]];
        end
        
        figure(1)
        set(gcf,'position',[50 50 1000 225])
        subplot(1,4,i)
        label = {sprintf('%i-%i',train)};
        label = [label {sprintf('%i-%i',test([1 2]))}];
        label = [label {sprintf('%i-%i',test([2 3]))}];
        label = [label {sprintf('%i-%i',test([1 3]))}];
        mkBow(toPlot,label)
        xlabel('Sequence similarity')
        ylabel('Rate map correlation (r)')
        hold on
        set(gca,'ylim',[-1 1])
        plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
    end
    saveFig(gcf,[root '/SimXSeq_MFRHists'],[{'tiff'} {'pdf'} {'jpeg'}])

    stl = mat2lag(across.simXseq(1:2,1:2,:));
%     ltl = mat2lag(across.simXseq(3:5,3:5,:));
    ltl = [permute(across.simXseq(3,4,:),[3 1 2]) ...
        permute(across.simXseq(4,5,:),[3 1 2]) ...
        permute(across.simXseq(3,5,:),[3 1 2])]; 
    isGood = ~isnan(stl(:,1)); % Only include cells with STL vals
    [blah b] = sort(stl(isGood,1),'ascend');  
    plotSimXLag(ltl(isGood,:),[root '/SimXSeqLag_SortedBy_ShortTermLag'],b)
    
    
    stl = mat2lag(across.simXseq(4:5,4:5,:));
%     ltl = mat2lag(across.simXseq(3:5,3:5,:));
    ltl = [permute(across.simXseq(1,2,:),[3 1 2]) ...
        permute(across.simXseq(2,3,:),[3 1 2]) ...
        permute(across.simXseq(1,3,:),[3 1 2])]; 
    isGood = ~isnan(stl(:,1)); % Only include cells with STL vals
    [blah b] = sort(stl(isGood,1),'ascend');  
    plotSimXLag(ltl(isGood,:),[root '/SimXSeqLag_SortedBy_ShortTermLag_Sort_4-5'],b)
    
    figure
    set(gcf,'position',[50 50 250 300])
    subplot(2,1,1)
    smfr = nanmean(across.mfr(isGood,1:12),2);
    gmfr = [{smfr(b(1:floor(length(b)./2)))} ...
        {smfr(b(floor(length(b)./2)+1:end))}];
    [blah matchedInds] = compHist(gmfr,[0:0.0025:0.08]);
    ylabel('Proportion')
    xlabel('Mean firing rate')
    subplot(2,1,2)
    compHist([{gmfr{1}(matchedInds(:,1))} ...
        {gmfr{2}(matchedInds(:,2))}],[0:0.0025:0.08]);
    ylabel('Proportion')
    xlabel('Mean firing rate')
    saveFig(gcf,[root '/SimXSeq_MFRHists'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    figure
    set(gcf,'position',[50 50  300 300])
    toPlot = ltl(isGood,:);
    gltl = [{toPlot(b(1:floor(length(b)./2)),:)} ...
        {toPlot(b(floor(length(b)./2)+1:end),:)}];
    toPlot = [{gltl{1}(matchedInds(:,1),1)} {gltl{2}(matchedInds(:,2),1)}; ...
        {gltl{1}(matchedInds(:,1),2)} {gltl{2}(matchedInds(:,2),2)}; ...
        {gltl{1}(matchedInds(:,1),3)} {gltl{2}(matchedInds(:,2),3)}];
    mkBow(toPlot')
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
    xlabel('Lag')
    ylabel('Similarity (r)')
    saveFig(gcf,[root '/SimXSeqLag_SortedBy_ShortTermLag_MatchedMFR'],[{'tiff'} {'pdf'} {'jpeg'}])
end