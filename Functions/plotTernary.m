function h = plotTernary(d,val,colors,labels,doBackground)
    
    if nargin < 5 || isempty(doBackground)
        doBackground = true;
    end

    if nargin < 2 || isempty(val)
        val = ones(length(d(:,1))).*20; 
    end
    
    if nargin < 3 || isempty(colors)
        colors = hsv(length(d(:,1)));
    end

%     nd = bsxfun(@rdivide,d,nansum(d,2))
    nd = d;
    angs = 0:360./length(d(1,:)):360-360./length(d(1,:));
    bangs = 0:360./(length(d(1,:)).*2):360-360./(length(d(1,:)).*2);
    
    x = nansum(bsxfun(@times,cosd(angs),nd),2);
    y = nansum(bsxfun(@times,sind(angs),nd),2);

    hold on
    
    if doBackground
        patch(sind(bangs([1:end 1])),cosd(bangs([1:end 1])),[0.5 0.5 0.5])
        plot(sind(bangs([1:end 1])),cosd(bangs([1:end 1])),...
            'color',[0.75 0.75 0.75],...
            'linestyle','-','linewidth',1.5)
        for i = [0.25:0.25:0.75]
            plot(i.*sind(bangs([1:end 1])),i.*cosd(bangs([1:end 1])),'color',[0.75 0.75 0.75],...
            'linestyle',':','linewidth',1)
        end
        plot([zeros(1,length(angs)); sind(angs)],[zeros(1,length(angs)); cosd(angs)], ...
            'color',[0.75 0.75 0.75],...
            'linestyle',':','linewidth',1)
    end
    
    h = scatter(y,x,val,colors,'filled','markeredgecolor','k');
    set(gca,'xlim',[-1.05 1.05]);
    axis equal
%     set(gca,'xlim',[-1 1],'ylim',[-1 1])
    axis off
    
    if nargin > 3 && ~isempty(labels)
        text(sind(angs).*1.15,cosd(angs).*1.15,labels,'horizontalalignment','center',...
            'verticalalignment','middle','fontname','arial','fontsize',11);
%         text(sind(angs+angs(2)./2).*1.15,cosd(angs+angs(2)./2).*1.15,labels,'horizontalalignment','center',...
%             'verticalalignment','middle','fontname','arial','fontsize',11);
    end
end