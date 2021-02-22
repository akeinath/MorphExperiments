function basicsOut(allR2Penalty,value,allTR2,tr2Thresh,root,outLabel)
    groupLabels = [{'Time'} {'Attractor'} {'Geometry'} {'Unclassified'}];
    toPlot = repmat({[]},[1 4]);
    [a b] = nanmax(allR2Penalty,[],2);
    b(isnan(allTR2)) = nan;
    for i = 1:3
        toPlot{i} = value(b==i & allTR2 > tr2Thresh);
    end
    toPlot{4} = value(allTR2 <= tr2Thresh);
    
    figure
    set(gcf,'position',[50 50 225 250],'color','w')
    mkBow(toPlot,groupLabels);
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    set(gca,'outerposition',[0 0 1 0.9])
    set(gca,'ylim',[nanmin(0,nanmin(get(gca,'ylim'))) nanmax(get(gca,'ylim'))])
    xlabel('Cell group')
    ylabel('Median within-session SHC (r)')
    saveFig(gcf,[root '/Overall/' outLabel],[{'pdf'} {'tiff'}])
    
    a = cellfun(@nanmean,toPlot);
    b = cellfun(@nanstd,toPlot)./sqrt(cellfun(@numel,toPlot));
    c = [a; b];
    fid = fopen(['Stats_' outLabel '.txt'],'w');
    fprintf(fid,'\n\tM+/-SEM: %0.4f +/- %0.6f',c);
    fprintf(fid,'\n\n\t\tSigned Rank vs 0\n');
    for i = 1:4
        [pval h stats] = signrank(toPlot{i},0);
        fprintf(fid,['\n(Group %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
            [i stats.signedrank stats.zval pval]);
    end
    fprintf(fid,'\n\n\t\tRank sum pairwise comparisons\n');
    for i = 1:4
        for j = i+1:4
            [pval h stats] = ranksum(toPlot{i},toPlot{j});
            fprintf(fid,['\n(' groupLabels{i} ' vs ' groupLabels{j} '): W = %.0f, Z = %0.2f, p = %.2e'], ...
                [stats.ranksum stats.zval pval]);
        end
    end
    fclose(fid);
end


















