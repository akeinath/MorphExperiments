function bv = bin(v,binStep)
    
    stack = [];
    for i = 1:binStep
        tmp = v(:,i:binStep:end,:);
        while ~isempty(stack) && length(tmp(1,:,1,1)) < length(stack(1,:,1,1))
            tmp = cat(2,tmp,zeros([size(tmp(:,1,:,1)) 1]));
        end
        stack = cat(4,stack,tmp);
    end
    bv = nanmean(stack,4);
    bv(:,end,:) = [];
end
