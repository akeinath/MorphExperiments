function mkBowSessions(doValue)
    set(gca,'ylim',[0 32])
    mkBow(doValue,[{'Time'} ...
        {'Attractor'} {'Geometry'}])
    ylabel('Tracked Sessions')
    xlabel('Cell group')
    apvals = nan(size(doValue));
    astats = nan(size(doValue));
    iter = 0;
    for i = 1:length(doValue)
        text(i,2,sprintf('n = %.0f',length(doValue{i})),...
                'horizontalalignment','center','color','k','fontsize',7)
        for j = i+1:length(doValue)
            iter = iter+1;
            [apvals(i,j) h stats] = ranksum(doValue{i},doValue{j});
            astats(i,j) = stats.ranksum;
            h = plot([i,j],ones(1,2).*(30+iter.*3),'color','k');
            if apvals(i,j) < 0.05
                set(h,'linewidth',2);
            end            
            text(nanmean([i j]),(30+iter.*3)+1.3,sprintf('W=%.0f, p=%.1e',[astats(i,j) apvals(i,j)]),...
                    'horizontalalignment','center','color','k','fontsize',7)
        end
    end
end