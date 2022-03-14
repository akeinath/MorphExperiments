function pairwiseNonparametrics(toPlot,root)

    
    if ~iscell(toPlot)
        tmp = repmat({[]},[1 length(toPlot(1,:))]);
        for i = 1:length(tmp(1,:))
            tmp{i} = toPlot(:,i);
        end
        toPlot = tmp;
    end

    a = cellfun(@nanmean,toPlot);
    b = cellfun(@nanstd,toPlot)./sqrt(cellfun(@numel,toPlot));
    c = [a; b; cellfun(@nanstd,toPlot)];
    checkP([root '_stats_nonparametric.txt']);
    fid = fopen([root '_stats_nonparametric.txt'],'w');
    fprintf(fid,'\n\tM+/-SEM +/-STD: %0.4f +/- %0.6f  +/- %0.6f',c);
    fprintf(fid,'\n');
    for i = 1:numel(toPlot)
        for j = i+1:numel(toPlot)
            if length(toPlot{i}) == length(toPlot{j})
                [pval h stats] = signrank(toPlot{i},toPlot{j});
                fprintf(fid,['\n(Signed rank %.0f vs %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
                    [i j stats.signedrank stats.zval pval]);
            end
            [pval h stats] = ranksum(toPlot{i},toPlot{j});
            fprintf(fid,['\n(Rank sum %.0f vs %.0f): W = %.0f, Z = %0.2f, p = %.2e'], ...
                [i j stats.ranksum stats.zval pval]);
            
        end
    end
    fclose(fid);
end