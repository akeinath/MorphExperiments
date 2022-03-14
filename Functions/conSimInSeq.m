function [seqSim simXlag] = conSimInSeq(sim,envs,envLabel)
    
    
    vals = sim.pearson;
    vals(31:32,:,:) = [];
    vals(:,31:32,:) = [];
    envs(31:32) = [];
    
    [a stereoseq] = ismember(envs,envLabel);
    [a b] = sort(stereoseq,'ascend');
    vals = vals(b,b,:);
    vals = nanmax(vals,permute(vals,[2 1 3]));
    
    isDiag = false(6,6);
    isDiag(1:7:numel(isDiag)) = true;
    isDiag = repmat(isDiag,size(vals(1,1,:)));
    
    numSeq = length(vals(:,1,1))./6;
    seqSim = nan(numSeq,numSeq,length(vals(1,1,:)));
    for i = 1:numSeq
        for j = i+1:numSeq
            tmp = vals(i:numSeq:end,j:numSeq:end,:);
            seqSim(i,j,:) = nanmean(reshape(tmp(isDiag),[6 length(tmp(1,1,:))]),1);
        end
    end
    
    simXlag = mat2lag(seqSim);
end