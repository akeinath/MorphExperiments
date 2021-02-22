function h = plotTernary(d,labels,doBackground)
    
    if nargin < 3 || isempty(doBackground)
        doBackground = true;
    end

    nd = d;
    angs = 0:360./length(d(1,:)):360-360./length(d(1,:));
    
    x = nansum(bsxfun(@times,cosd(angs),nd),2);
    y = nansum(bsxfun(@times,sind(angs),nd),2);

    hold on
    
    ps = 0.01;
    
    bx = floor((x+1)./ps)+1;
    by = floor((y+1)./ps)+1;
    m = zeros(floor(2./ps)+1);
    for i = 1:nanmax(bx)
        for j = 1:nanmax(by)
            m(i,j) = nansum(bx==i&by==j);
        end
    end
    
    m = imfilter(m,fspecial('gauss',[15 15],4),'same');
    imagesc(-1:ps:1,-1:ps:1,m./nansum(m(:)));
    colormap(hot)
    set(gca,'ydir','normal');
    
    if doBackground
        plot(sind(angs([1:end 1])),cosd(angs([1:end 1])),...
            'color',[0.75 0.75 0.75],...
            'linestyle','-','linewidth',1.5)
        for i = [0.25:0.25:0.75]
            plot(i.*sind(angs([1:end 1])),i.*cosd(angs([1:end 1])),'color',[0.75 0.75 0.75],...
            'linestyle',':','linewidth',1)
        end
        plot([zeros(1,length(angs)); sind(angs)],[zeros(1,length(angs)); cosd(angs)], ...
            'color',[0.75 0.75 0.75],...
            'linestyle',':','linewidth',1)
    end
    
    set(gca,'xlim',[-1.1 1.1],'color','k');
    axis equal
%     set(gca,'xlim',[-1 1],'ylim',[-1 1])
    axis off
    
    
    if doBackground
        c = colorbar('south','color','w','AxisLocation','out');
        c.Label.String = 'Proportion of cells';
    end
    
    if nargin > 1 && ~isempty(labels)
        text(sind(angs).*1.15,cosd(angs).*1.15,labels,'horizontalalignment','center',...
            'verticalalignment','middle','fontname','arial','fontsize',11,'color','w');
    end
end