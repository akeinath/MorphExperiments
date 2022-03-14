function [r2penalty tr2 shuffle_r2penalty shuffle_tr2 pop_r2penalty pop_tr2] = classFactorLoading(vals,RDMs,minimum_pairwise_comparisons,root)

    [r2penalty tr2] = classNFactor(cat(3,vals,nanmean(vals,3)),RDMs, ...
        minimum_pairwise_comparisons);
    
    pop_tr2 = tr2(end);
    pop_r2penalty = r2penalty(end,:);
    tr2(end) = [];
    r2penalty(end,:) = [];

    % do Shuffled;
    tmp = vals;
    goodCells = nansum(~isnan(reshape(tmp,[numel(tmp(:,:,1)) length(tmp(1,1,:))]))) >= ...
        minimum_pairwise_comparisons;
    tmp(:,:,~goodCells) = nan;
    for i = 1:length(tmp(:,1,1))
        for j = 1:length(tmp(1,:,1))
            isGood = find(~isnan(tmp(i,j,:)));
            pig = isGood(randperm(length(isGood)));
            tmp(i,j,isGood) = tmp(i,j,pig);
        end            
    end

    [shuffle_r2penalty shuffle_tr2] = classNFactor(tmp,RDMs, ...
        minimum_pairwise_comparisons);
end