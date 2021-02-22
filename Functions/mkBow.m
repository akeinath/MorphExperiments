function h = mkBow(d,xl,varargin)
    if nargin<2 || isempty(xl)
        xl = 1:length(d(1,:));
    end

    if ~iscell(d)
        newD = repmat({[]},[1 length(d(1,:))]);
        for i = 1:length(d(1,:))
            newD{i} = d(:,i);
        end
        d = newD;
    end
    
    dl = cellfun(@size,d,repmat({2},size(d)));
    if dl>1
        d = d';
        new = repmat({[]},[length(d) nanmax(dl)]);
        for j = 1:length(d)
            for k = 1:nanmax(dl)
                try
                    new{j,k} = d{j}(:,k);
                end
            end
        end
        d = new;
    end
    
    groups = length(d(:,1));
    ticks = length(d(1,:));
    set(gca,'xlim',[0.5 ticks+0.5],'xtick',[1:ticks],'xticklabel',xl,'fontname','arial',...
        'fontsize',10,'fontweight','normal')
    hold on
    
    w = 0.6;
    wpg = w./groups;
    groupColor = [1 0.6 0.2; 0.2 0.6 1];
    dip = 0.025;

    for i = 1:groups
        for j = 1:ticks
            if isempty(d{i,j})
                continue
            end
            sVals = sort(d{i,j});
            sVals(isnan(sVals)) = [];
            
            cutoff = [nanmean(sVals)-nanstd(sVals) nanmean(sVals)+nanstd(sVals)];
            cutoff = [-inf inf];

            lqr = sVals(floor(length(sVals).*0.25));
            uqr = sVals(floor(length(sVals).*0.75));
            mqr = sVals(floor(length(sVals).*0.5));
            
            doC = groupColor(mod((j-1),length(groupColor(:,1)))+1,:);
            if groups >1
                shift = -0.2+(i-1).*0.4;
            else
                shift = 0;
            end
            plot(j.*ones(2,2)+shift,[max(cutoff(1),sVals(1)) lqr; min(cutoff(2),sVals(end)) uqr]', ...
                'color',doC,'linewidth',1.5,'linestyle','-');
            h(i,j) = patch(j+(w./2./groups).*[-1 1 1 0.5 -0.5 -1]+shift,...
                [lqr lqr mqr-dip.*range(get(gca,'ylim')) ...
                mqr mqr mqr-dip.*range(get(gca,'ylim'))], ...
                'w','facecolor',doC.*0.5+0.5,'edgecolor',doC,'linewidth',1.5);
            h(i,j) = patch(j+(w./2./groups).*[-1 1 1 0.5 -0.5 -1]+shift,...
                [uqr uqr mqr+dip.*range(get(gca,'ylim')) ...
                mqr mqr mqr+dip.*range(get(gca,'ylim'))],...
                'w','facecolor',doC.*0.5+0.5,'edgecolor',doC,'linewidth',1.5); %doC./2+0.5
            
%             plot(mean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),...
%                 sVals([sVals<cutoff(1)|sVals>cutoff(2)]),'linestyle','none',...
%                 'marker','d','markersize',2.5,'markerfacecolor','k','color','k','markeredgecolor','k')
            
%             h(i,j) = patch([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))],...
%                 [sVals(floor(length(sVals).*0.25)) sVals(floor(length(sVals).*0.25)) ...
%                 sVals(ceil(length(sVals).*0.75)) sVals(ceil(length(sVals).*0.75))],...
%                 groupColor(ceil(i),:),'edgecolor',groupColor(ceil(i),:),'linewidth',2);
%             
%             plot(([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),...
%                 [nanmedian(sVals) nanmedian(sVals)],'color','k','linewidth',2);
        end
    end
end