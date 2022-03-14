function h = mkAcrossLine(d,xl,groupColor)    

    % mkLine([{randn(400,10)}; {randn(400,10)}]) <-------- example

    if ~iscell(d)
        d = {d};
%         tmp = repmat({[]},[1 length(d(1,:))]);
%         for i = 1:length(d(:,1))
%             tmp{i} = d(:,i);
%         end
%         d = tmp;
    end
    
     if nargin<2 || isempty(xl)
        xl = 1:length(d(1,:));
    end
    
    set(gca,'xlim',[xl(1)-0.5 xl(end)+0.5],'xtick',xl,'xticklabel',xl,'fontname','arial',...
        'fontsize',9,'fontweight','normal')
    hold on
    
    w = 0.8;

    if nargin < 3 || isempty(groupColor)
        groupColor = [1 0.4 0.4; 0.4 0.4 1; 0.2 1 0.2; 0.75 0.75 0.75; 0.4 0.9 1];
    end

    for i = 1:length(d(:,1))
        doSpot = xl;
        
        m = cellfun(@nanmean,d(i,:));
        sd = cellfun(@nanstd,d(i,:));
        n = cellfun(@nansum,cellfun(@not,cellfun(@isnan,d(i,:),'uniformoutput',false),'uniformoutput',false));
        sem = sd./sqrt(n);
        
        quarts = [m-sem; m; m+sem];
        
        patch([doSpot(1:length(doSpot)-1); doSpot(1:length(doSpot)-1); ... 
            doSpot(2:length(doSpot)); doSpot(2:length(doSpot))],...
            [quarts(3,1:end-1); quarts(1,1:end-1); ...
            quarts(1,2:end); quarts(3,2:end)],groupColor(i,:),'edgecolor','none',...
            'facecolor',groupColor(i,:),'edgealpha',0,'facealpha',0.75)
        h(i) = plot(xl,quarts(2,:),'color',groupColor(i,:),'linestyle','-','linewidth',1.25);       
    end
    
%     if nanmin(get(gca,'ylim'))<0
%         plot(get(gca,'xlim'),[0 0],'linestyle','--','linewidth',1,'color','k');
%     end
end