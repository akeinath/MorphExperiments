function accuracies = decodeContext_SVM(gT,envs,stepSize,setSize)

    if nargin < 3 || isempty(stepSize)
        stepSize = 6;
    end
    if nargin < 4 || isempty(setSize)
        setSize = 6;
    end
    
    binSize = 60;
    
    cGT = [];
    for i = 1:length(gT)
        cGT = [cGT {chunkVec(gT{i},binSize)}];
    end
    
    allEnvAcc = [];
    allCGAcc = [];
    allSims = [];
    trainSet = bsxfun(@plus,[0:(setSize-1)],[1:stepSize:length(gT)-(setSize-1)]');
    for ti = 1:length(trainSet(:,1))
        test = setxor([1:length(gT)],trainSet(ti,:));
        
        isGroup = strcmp(cellfun(@(x) x(1),envs,'un',0),{'g'})';
        
        trainEnvs = envs(trainSet(ti,:));
        isTrainGroup = isGroup(trainSet(ti,:));
        atrVec = [];
        labels = [];
        for tri = 1:length(trainSet(ti,:))
            tmp = cGT{trainSet(ti,tri)};
            atrVec = cat(2,atrVec,tmp);
            labels = cat(2,labels,ones(1,length(tmp(1,:))).*isTrainGroup(tri));
        end
        
        allAcc = nan(1,length(test));
        for testi = 1:length(test)
            
            include = ~any(isnan(atrVec),2) & ~any(isnan(cGT{test(testi)}),2);

            outSVM = fitcsvm(atrVec(include,:)',labels(:));
            pred = predict(outSVM,cGT{test(testi)}(include,:)');
            allAcc(testi) = nanmean(pred==isGroup(test(testi)));
        end  

        lags = nanmin(abs(bsxfun(@minus,trainSet(ti,:)',test)),[],1);
        allCGAcc = [allCGAcc [lags; allAcc]];
    end
    
    blah = [1:6:nanmax(allCGAcc(1,:))];
    blah(end) = inf;
    lagGroup = [blah(1:end-1)' blah(2:end)'-1];
    accuracies = nan(length(lagGroup(:,1)),3);
    for i = 1:length(lagGroup(:,1))
        accuracies(i,2) = nanmean(allCGAcc(2,allCGAcc(1,:)>=lagGroup(i,1) & ...
            allCGAcc(1,:)<=lagGroup(i,2)));
    end
end