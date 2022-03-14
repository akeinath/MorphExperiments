function [conSim simXlag] = conComparisonInSeq(sim,envs,envLabel)

    vals = sim.pearson;
%     vals(repmat(~inc,[1 1 length(sim.pearson(1,1,:))])) = nan;
    
    seqLag = abs(bsxfun(@minus,ceil([1:length(vals(:,1,1))]./6), ...
        ceil([1:length(vals(:,1,1))]'./6)));
    vals(repmat(seqLag~=0,[1 1 length(vals(1,1,:))])) = nan;
    
    [a envType] = ismember(envs,envLabel);
    
    conSim = nan(length(envLabel),length(envLabel),length(sim.pearson(1,1,:)));
    for i = 1:length(envLabel)
        for j = i+1:length(envLabel)
            isGood = bsxfun(@times,envType==i,envType'==j) | ...
                bsxfun(@times,envType==j,envType'==i);
            tv = vals(repmat(isGood,[1 1 length(vals(1,1,:))]));
            tv = reshape(tv,[nansum(isGood(:)) length(vals(1,1,:))]);
            conSim(i,j,:) = nanmean(tv,1);
        end
    end
    
    simXlag = mat2lag(conSim);
end