function [c total_comps] = cellwiseMatCorr(a,b)

    a = reshape(a,[numel(a(:,:,1)) length(a(1,1,:))]);
    b = reshape(b,[numel(b(:,:,1)) length(b(1,1,:))]);

    c = nan(length(a(1,:)));
    total_comps = zeros(length(a(1,:)));
    for i = 1:length(a(1,:))
        for j = i; %1:length(a(1,:))
            isGood = ~isnan(a(:,i))&~isnan(b(:,j));
            total_comps(i,j) = sum(isGood);
            if sum(isGood)>1
                c(i,j) = corr(a(isGood,i),b(isGood,j));
            end
        end
    end
    c = c(1:length(c)+1:end);
    total_comps = total_comps(1:length(total_comps)+1:end);
end