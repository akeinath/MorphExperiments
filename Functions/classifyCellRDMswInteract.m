function [outR2Penalty tr2] = classifyCellRDMswInteract(sim,RDMs,root,doShuffle)
    
    minSessions = 28;

    if nargin < 4 || isempty(doShuffle)
        doShuffle = false;
    end
   
    rsim = reshape(sim,[numel(sim(:,:,1)) length(sim(1,1,:))]);
    rRDMs = reshape(RDMs,[numel(RDMs(:,:,1)) length(RDMs(1,1,:))]);

    rsim = atanh(rsim); % normalize
    
    % z-score for shuffling so baseline differences aren't messy
%     rsim = bsxfun(@rdivide,bsxfun(@minus,rsim,nanmean(rsim,1)),nanstd(rsim,[],1));
    
    if doShuffle
        isGood = (nansum(~isnan(rsim),1)>=minSessions);
%         isGood = true(1,length(rsim(1,:)));
        for i = 1:length(rsim(:,1))
            gi = find(isGood&~isnan(rsim(i,:)));
            if ~isempty(gi)
                sgi = gi(randperm(length(gi)));
                rsim(i,gi) = rsim(i,sgi);
            end
        end
    end
    
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
    
    threeWay = prod(rRDMs,2);
    
    tr2 = nan(length(rsim(1,:)),1);
    dropr2 = nan(length(rsim(1,:)),length(rRDMs(1,:)));
    for k = 1:length(rsim(1,:))
        if nansum(~isnan(rsim(:,k))) < minSessions
            continue
        end
        
        x = rsim(~isnan(rsim(:,k)),k);
        [B dev stats] = glmfit([rRDMs(~isnan(rsim(:,k)),:) ...
            twoWay(~isnan(rsim(:,k)),:) threeWay(~isnan(rsim(:,k)),:)],x,'normal');
        tr2(k) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        
        for j = 1:length(rRDMs(1,:))
            [B dev stats] = glmfit([rRDMs(~isnan(rsim(:,k)),[1:j-1 j+1:end]) ...
                twoWay(~isnan(rsim(:,k)),~any(usesRDMs==j,1))],x,'normal');
            dropr2(k,j) = 1 - nansum(stats.resid.^2)./nansum((x-nanmean(x)).^2);
        end
    end
    
    r2penalty = 1-bsxfun(@rdivide,dropr2,tr2);
    outR2Penalty = r2penalty;
    r2penalty  = r2penalty(~isnan(r2penalty(:,1)),:);
    
    % drop overall before returning
    outR2Penalty(1,:) = [];
    tr2(1,:) = [];
end