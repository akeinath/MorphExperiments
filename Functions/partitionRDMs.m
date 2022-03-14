function conSim = partitionRDMs(vals,envs,envLabel)

    if nargin < 3 || isempty(envLabel)
        envLabel = [{'g1'} {'g2'} {'g3'} {'sq3'} {'sq2'} {'sq1'}];
    end

    seqNum = ceil([1:length(vals(:,1,1))]./6);
    [a envType] = ismember(envs,envLabel);
    
    conSim = nan(length(envLabel),length(envLabel), ...
        nanmax(seqNum)-1,nanmax(seqNum)-1,length(vals(1,1,:)));
    for i = 1:length(envLabel)
        for j = i:length(envLabel)
            for k = 1:nanmax(seqNum)-1
                for q = k:1:nanmax(seqNum)-1
                    isGood = [bsxfun(@times,envType==i,envType'==j) | ...
                        bsxfun(@times,envType==j,envType'==i)] & ...
                        [bsxfun(@times,seqNum==k,seqNum'==q) | ...
                        bsxfun(@times,seqNum==q,seqNum'==k)];
                    tv = vals(repmat(isGood,[1 1 length(vals(1,1,:))]));
                    tv = reshape(tv,[nansum(isGood(:)) length(vals(1,1,:))]);
                    if nanmax(nansum(~isnan(tv),1))==1
                        conSim(i,j,k,q,:) = nanmean(tv,1);
                    end
                end
            end
        end
    end
end