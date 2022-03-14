function plotSimXLag(sl,root,sortCrit,numDivs)


    if nargin < 3 || isempty(sortCrit)
        [a sortCrit] = sort(sl(:,1),'descend');
    end
    ssl = sl(sortCrit,:);

    if nargin < 4 || isempty(numDivs)
        numDivs = 2;
    end

    
    bounds = 1:ceil(length(ssl(:,1))./numDivs):length(ssl(:,1));
    bounds(end+1) = length(ssl(:,1))+1;
    gd = repmat({[]},[1 length(bounds)-1]);
    for i = 1:length(gd)
        gd{i} = ssl(bounds(i):bounds(i+1)-1,:);
    end
    
    %%% Make figure
    
    figure
    set(gcf,'position',[50 50  300 300])
%     h = mkLine(gd',[1:length(gd{1}(1,:))],transcm(length(gd)));
    
    h = mkBow(gd',[1:length(gd{1}(1,:))],transcm(length(gd)));
%     set(gca,'xlim',[0.75 length(gd{1}(1,:)) + 0.25], ...
%         'ylim',[-0.1 0.8],'fontname','arial','fontweight','normal')
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
%     plot([1.5 1.5],get(gca,'ylim'),'color','k','linestyle','--')
%     legend(h,[{'0-20%'} {'20-40%'} {'40-60%'} {'60-80%'} {'80-100%'}],'location','northeast')
    
%     tmp = [0:1./numDivs:1];
%     tmp = round([tmp(1:end-1)' tmp(2:end)'].*100);
%     tl = [];
%     for i = 1:length(tmp(:,1))
%         tl = [tl {[num2str(tmp(i,1)) '-' num2str(tmp(i,2)) '%']}];
%     end
%     set(gca,'xticklabel',tl);
    xlabel('Lag')
    ylabel('Similarity (r)')
    saveFig(gcf,[root],[{'tiff'} {'pdf'} {'jpeg'}])
end