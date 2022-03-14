function [simXlag valsXlag] = mat2lag(seqSim)
    
    sseqSim = seqSim;
    
    mask = false(size(seqSim(:,:,1)));
    mask(triu(true(size(mask)),1)) = true;
    for lag = 1:length(seqSim(1,:,1))-1
        sseqSim(1+lag,:,:) = circshift(seqSim(1+lag,:,:),[0 -lag 0]);
        mask(1+lag,:) = circshift(mask(1+lag,:,:),[0 -lag 0]);
    end
    sseqSim(~mask) = nan;
    tmp = permute(nanmean(sseqSim,1),[3 2 1]);
    simXlag = tmp(:,2:end);
    
    valsXlag = reshape(sseqSim(repmat(mask,[1 1 length(sseqSim(1,1,:))])), ...
        [nansum(mask(:)) length(sseqSim(1,1,:))])';
end