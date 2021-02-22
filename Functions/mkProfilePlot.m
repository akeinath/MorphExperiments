function mkProfilePlot(d)
    
    colors = [0.9 0.2 0.2; 0.2 0.2 0.9; 0.2 0.7 0.2];
    for k = 1:length(d)        
        m = nanmean(d{k});
        se = nanstd(d{k})./sqrt(nansum(~isnan(d{k})));
        x = [0:31];
        hold on
        plot(x,m+se,'color',colors(k,:).*0.5+0.5);
        plot(x,m-se,'color',colors(k,:).*0.5+0.5);
        plot(x,m,'color',colors(k,:));
        hold on
    end
end


















