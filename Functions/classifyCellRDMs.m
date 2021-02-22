function classifyCellRDMs(sim,RDMs,root)
    
   
    rsim = reshape(sim,[numel(sim(:,:,1)) length(sim(1,1,:))]);
    rRDMs = reshape(RDMs,[numel(RDMs(:,:,1)) length(RDMs(1,1,:))]);
    
    rsim = [nanmedian(rsim,2) rsim];
    
    % normalize RDMs
    for j = 1:length(rRDMs(1,:))
        rRDMs(:,j) = (rRDMs(:,j)-nanmin(rRDMs(:,j)));
        rRDMs(:,j) = (rRDMs(:,j)./nanmax(rRDMs(:,j)));
    end
    
    tBs = nan(length(rsim(1,:)),length(rRDMs(1,:))+1);
    tr2 = nan(length(rsim(1,:)),1);
    dropr2 = nan(length(rsim(1,:)),length(rRDMs(1,:)));
    for k = 1:length(rsim(1,:))
        if nansum(~isnan(rsim(:,k))) < 50
            continue
        end
        
        x = rsim(~isnan(rsim(:,k)),k);
        [B dev stats] = glmfit(rRDMs(~isnan(rsim(:,k)),:),x,'normal');
        tBs(k,:) = B;
        tr2(k) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        
        for j = 1:length(rRDMs(1,:))
            [B dev stats] = glmfit(rRDMs(~isnan(rsim(:,k)),[1:j-1 j+1:end]),x,'normal');
            dropr2(k,j) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        end
    end
%     figure
%     set(gcf,'position',[50 50 75.*length(RDMs(1,1,:)) 225])
%     mkWhisker(tBs(:,2:end));
%     hold on
%     plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--');
    
    r2penalty = 1-bsxfun(@rdivide,dropr2,tr2);
    r2penalty  = r2penalty(~isnan(r2penalty(:,1)),:);
    
    plotWeightClass(r2penalty,tr2,[root '/Classification/Dropout_r2.gif'],true)
end