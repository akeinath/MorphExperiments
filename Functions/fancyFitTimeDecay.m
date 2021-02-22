function [ap val] = fancyFitTimeDecay(sim)

%     sim = atanh(sim);

    rsim = [nan(size(sim)) sim];
    for i = 1:length(sim(:,:,1))
        rsim(i,:,:) = circshift(rsim(i,:,:),[0 -(i-1)]);
    end
    val = permute(nanmean(rsim(:,33:end,:),1),[3 2 1]);
    
    ap = nan(length(rsim(1,1,:)),2);
    x = repmat(abs([-31:32]),[32 1]);
    for k = 1:length(rsim(1,1,:))
        y = rsim(:,:,k);
        [p S] = polyfit(x(~isnan(y)),y(~isnan(y)),1);
        ap(k,:) = fliplr(p);
        
    end
end