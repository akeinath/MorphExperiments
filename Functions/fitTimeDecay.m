function [ap val] = fitTimeDecay(sim)

%     sim = atanh(sim);

    rsim = [nan(size(sim)) sim];
    for i = 1:length(sim(:,:,1))
        rsim(i,:,:) = circshift(rsim(i,:,:),[0 -(i-1)]);
    end
    rsim(:,1:length(sim(1,:,1)),:) = [];
    val = permute(nanmean(rsim,1),[3 2 1]);
    
    ap = nan(length(val(:,1)),2);
    x = [1:length(val(1,:))]-1;
    for k = 1:length(val(:,1))
        y = val(k,:);
        [p S] = polyfit(x(~isnan(y)),y(~isnan(y)),1);
        ap(k,:) = fliplr(p);
    end
end