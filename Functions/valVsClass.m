function atab = valVsClass(vals,r2penalty,tr2,labels)

    if nargin < 4 || isempty(labels)
        labels = repmat({'Blah'},[1 length(r2penalty(1,:))]);
    end

    exclude = any(isnan(tr2),2);
    vals(exclude,:) = [];
    r2penalty(exclude,:) = [];
    tr2(exclude,:) = [];

    colors = bsxfun(@times,[1 0 0],r2penalty(:,1));
    colors = colors+bsxfun(@times,[0 0 1],r2penalty(:,2));
    colors = colors+bsxfun(@times,[0 1 0],r2penalty(:,3));
    colors(colors<0) = 0;
    colors = colors.^(1./4);
    colors = [colors./nanmax(colors(:))];

    %%%% Correlation with firing properties

    gd = repmat({[]},[2 length(r2penalty(1,:))]);
    
    figure
    set(gcf,'position',[50 50 300.*(length(r2penalty(1,:))+1) 300])
    for j = 1:length(r2penalty(1,:))
        include = (r2penalty(:,j)==nanmax(r2penalty,[],2));
        subplot(1,length(r2penalty(1,:))+1,j)
        scatter(vals(include),tr2(include), ...
            30,colors(include,:),'filled')
        h = lsline;
        set(h,'color',[0 0.75 0.75],'linewidth',2)
        title(labels{j});
        ylabel('Full model r^2')
        xlabel('Value')
        set(gca,'ylim',[0 1]);
        
        set(gca,'xlim',[-0.5 1]);
        axis square
        
        gd{1,j} = vals(include);
        gd{2,j} = tr2(include);
        
        [rval pval] = corr(vals(include),tr2(include));
        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.9;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.1;
        text(textX,textY,sprintf(['r = %.3f\np = %.2e\nn = %.0f'],[rval pval nansum(include)]), ...
            'fontweight','normal','fontname','arial','fontsize',10);
    end
    subplot(1,length(r2penalty(1,:))+1,length(r2penalty(1,:))+1)
    scatter(vals,tr2,30,colors,'filled')
    h = lsline;
    set(h,'color',[0 0.75 0.75],'linewidth',2)
    title('Overall')
    ylabel('Full model r^2')
    xlabel('Value')
    set(gca,'ylim',[0 1]);
    axis square
    
    [rval pval] = corr(vals,tr2);
    textX = get(gca,'xlim');
    textX = textX(2) - [textX(2)-textX(1)].*0.9;
    textY = get(gca,'ylim');
    textY = textY(2) - [textY(2)-textY(1)].*0.1;
    text(textX,textY,sprintf(['r = %.3f\np = %.2e\nn = %.0f'],[rval pval length(tr2)]), ...
        'fontweight','normal','fontname','arial','fontsize',10);
    
    [h atab ctab stats] = aoctool(cat(1,gd{1,:}),cat(1,gd{2,:}),[ones(length(gd{1,1}),1).*1; ...
        ones(length(gd{1,2}),1).*2; ones(length(gd{1,3}),1).*3],[],[],[],[],'off','separate lines');
end