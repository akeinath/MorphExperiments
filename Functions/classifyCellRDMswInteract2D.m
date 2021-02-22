function classifyCellRDMswInteract2D(sim,RDMs,root)
    
   
    rsim = reshape(sim,[numel(sim(:,:,1)) length(sim(1,1,:))]);
    rRDMs = reshape(RDMs,[numel(RDMs(:,:,1)) length(RDMs(1,1,:))]);
    
    rsim = [nanmedian(rsim,2) rsim];
    
    % normalize RDMs
    for j = 1:length(rRDMs(1,:))
        rRDMs(:,j) = (rRDMs(:,j)-nanmin(rRDMs(:,j)));
        rRDMs(:,j) = (rRDMs(:,j)./nanmax(rRDMs(:,j)));
    end
    
    usesRDMs = [];
    twoWay = [];
    for j = 1:length(rRDMs(1,:))
        for k = j+1:length(rRDMs(1,:))
            twoWay = [twoWay rRDMs(:,j).*rRDMs(:,k)];
            usesRDMs = [usesRDMs [j; k]];
        end
    end
    
    tr2 = nan(length(rsim(1,:)),1);
    dropr2 = nan(length(rsim(1,:)),length(rRDMs(1,:)));
    for k = 1:length(rsim(1,:))
        if nansum(~isnan(rsim(:,k))) < 50
            continue
        end
        
        x = rsim(~isnan(rsim(:,k)),k);
        [B dev stats] = glmfit([rRDMs(~isnan(rsim(:,k)),:) ...
            twoWay(~isnan(rsim(:,k)),:)],x,'normal'); % twoWay(~isnan(rsim(:,k)),:)
        tr2(k) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        
        for j = 1:length(rRDMs(1,:))
            [B dev stats] = glmfit([rRDMs(~isnan(rsim(:,k)),[1:j-1 j+1:end]) ...
                twoWay(~isnan(rsim(:,k)),~any(usesRDMs==j,1))],x,'normal');
            dropr2(k,j) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        end
    end
    
    r2penalty = 1-bsxfun(@rdivide,dropr2,tr2);
    r2penalty  = r2penalty(~isnan(r2penalty(:,1)),:);
    
    
    colors = bsxfun(@times,[0 0.25 1],r2penalty(:,2));
    colors = colors+bsxfun(@times,[1 0 0],r2penalty(:,1));
    colors(colors<0) = 0;
    colors = colors.^(1./4);
    colors = [colors./nanmax(colors(:))];

    overallR2Penalty = r2penalty(1,:);
    overallTR2 = tr2(1,:);
    overallColor = colors(1,:);
    r2penalty(1,:) = [];
    tr2(1,:) = [];
    colors(1,:) = [];

    [blah best] = nanmax(tr2(~isnan(tr2)));
    
    figure
    set(gcf,'position',[50 50 400 200])
    hold on
    scatter(r2penalty(:,1)-r2penalty(:,2),0.10.*randn(length(r2penalty),1), ...
        (200.*tr2(~isnan(tr2))),colors,'filled','markeredgecolor','w');
    set(gca,'ylim',[-0.5 1],'xlim',[-1.1 1.1])
    scatter(overallR2Penalty(:,1)-overallR2Penalty(:,2),0, ...
        (500.*overallTR2),overallColor,'filled','marker','p','markeredgecolor','k');
    text(0,0.9,['Population r^2 = ' num2str(overallTR2)],...
        'horizontalalignment','center','color',[0.25 0.25 0.25],'fontsize',11,'fontweight','bold')  
    text(0,0.6,['Max Cell r^2 = ' num2str(nanmax(tr2))],...
        'horizontalalignment','center','color',[0.25 0.25 0.25],'fontsize',11,'fontweight','bold')  
    xlabel('Time - Attractor (dropout \Deltar^2%)');
    saveFig(gcf,[root '/Classification/Dropout_withInteractions_2D_r2'],[{'pdf'} {'tiff'}])
end