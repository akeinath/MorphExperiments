function addCorr(x,y,relX,relY)
    
    if nargin < 3 || isempty(relX)
        relX = 0.05;
    end
    if nargin < 4 || isempty(relY)
        relY = 0.95;
    end

    [rval pval] = corr(x,y);
    textX = get(gca,'xlim');
    textX = textX(1) + [textX(2)-textX(1)].*relX;
    textY = get(gca,'ylim');
    textY = textY(1) + [textY(2)-textY(1)].*relY;
    text(textX,textY,sprintf(['r = %.3f, p = %.2e'],[rval pval]), ...
        'fontweight','normal','fontname','arial','fontsize',9);
end