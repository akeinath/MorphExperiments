function lagVals = getCrossvalStab(ss)
    lagVals = repmat({[]},[4 4]);
    lag = abs(bsxfun(@minus,1:5,[1:5]'));
    for i = 1:length(ss(1,:,1))-1
        for j = i+1:length(ss(1,:,1))
            trainLag = lag(i,j);
            
            sts = permute(ss(i,j,:),[3 1 2]);
            
            doSim = setdiff([1:5],[i j]);
            lts = ss(doSim,doSim,:);
            isGood = ~all(isnan(lts),3);
            risGood = repmat(isGood,[1 1 length(lts(1,1,:))]);
            lts = permute(reshape(lts(risGood),[nansum(isGood(:)) length(lts(1,1,:))]),[2 1]);
            tlag = lag(doSim,doSim);
            tlag = tlag(isGood)';
            
            for k = 1:length(lts(1,:))
                tmp = [sts(:) lts(:,k)];
                goodVals = ~any(isnan(tmp),2);
                lagVals{trainLag,tlag(k)} = [lagVals{trainLag,tlag(k)}; tmp(goodVals,:)];
            end
        end
    end
end