function lagVals = getConVsStab(conMAE,ss)

    lagVals = repmat({[]},[1 4]);
    for i = 1:length(conMAE(1,:))
        doSim = [1:i-1 i+1:length(conMAE(1,:))];
        ts = ss(doSim,doSim,:);
        lag = abs(bsxfun(@minus,doSim,doSim'));
        
        tcm = repmat(permute(conMAE(:,i),[3 2 1]),[length(doSim) length(doSim)]);
        tlm = repmat(lag(:),[1 1 length(tcm(1,1,:))]);
        vals = [tcm(:) ts(:) tlm(:)];
        
        goodVals = vals(~any(isnan(vals),2),:);
        
        for j = 1:4
            lagVals{j} = [lagVals{j}; goodVals(j==round(goodVals(:,3)),1:2)];
        end
    end
end