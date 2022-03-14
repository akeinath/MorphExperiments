function out = matchNumRegisteredCells(in)
    counts = nansum(~isnan(in),3);
    isGood = triu(true(size(counts)),1);
    matchNum = nanmin(counts(isGood));
    
    out = nan(size(in));
    for si = 1:length(in(:,:,1))
        for sj = si+1:length(in(:,:,1))
            goodCells = find(~isnan(in(si,sj,:)));
            goodCells = goodCells(randperm(length(goodCells)));
            out(si,sj,goodCells(1:matchNum)) = in(si,sj,goodCells(1:matchNum));
        end
    end
end