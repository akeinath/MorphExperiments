function [Bs tr2] = classNFactor(sim,RDMs,minSessions)
   
    if ~iscell(RDMs)
        tmp = repmat({[]},[1 length(RDMs(1,1,:))]);
        for i = 1:length(RDMs(1,1,:))
            tmp{i} = RDMs(:,:,i);
        end
        RDMs = tmp;
    end
    
    rsim = reshape(sim,[numel(sim(:,:,1)) length(sim(1,1,:))]);
    
    rRDMs = repmat({[]},[1 numel(RDMs)]);
    for i = 1:numel(RDMs)
        rRDMs{i} = reshape(RDMs{i},[numel(RDMs{i}(:,:,1)) length(RDMs{i}(1,1,:))]);
    end
    
    rsim = atanh(rsim); % normalize
    
    % z-score for shuffling so baseline differences aren't messy
%     rsim = bsxfun(@rdivide,bsxfun(@minus,rsim,nanmean(rsim,1)),nanstd(rsim,[],1));
    
    % normalize RDMs
    for i = 1:numel(rRDMs)
        for j = 1:length(rRDMs{i}(1,:))
            rRDMs{i}(:,j) = (rRDMs{i}(:,j)-nanmin(rRDMs{i}(:,j)));
            rRDMs{i}(:,j) = (rRDMs{i}(:,j)./nanmax(rRDMs{i}(:,j)));
        end
    end
    
    % Make interaction RDMs, determine which groups, and concatonatve
    grdms = cat(2,rRDMs{:});
    aRDMs = grdms;
    gid = cellfun(@size,rRDMs,repmat({2},[1 numel(rRDMs)]));
    isGroup = [];
    for i = 1:length(gid)
        if i == 1
            isGroup(1:gid(i)) = i;
        else
            isGroup = [isGroup ones(1,gid(i)).*i];
        end
    end

    tmp = false(numel(rRDMs),length(isGroup));
    for i = 1:length(isGroup)
        tmp(isGroup(i),i) = true;
    end
    gid = tmp;
    
    tr2 = nan(length(rsim(1,:)),1); % full model r-squared
    Bs = nan(length(rsim(1,:)),length(aRDMs(1,:))); % dropout r-squared
    for k = 1:length(rsim(1,:))
        if nansum(~isnan(rsim(:,k))) < minSessions
            continue
        end
        
        isGood = ~isnan(rsim(:,k));
        x = rsim(isGood,k);
        [B dev stats] = glmfit([aRDMs(isGood,:)],x,'normal');
        tr2(k) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        Bs(k,:) = B(2:end)';
    end
end