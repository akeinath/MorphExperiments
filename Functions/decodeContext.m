function accuracies = decodeContext(sim,envs,stepSize,setSize)

    if nargin < 3 || isempty(stepSize)
        stepSize = 6;
    end
    if nargin < 4 || isempty(setSize)
        setSize = 6;
    end
    
    allEnvAcc = [];
    allCGAcc = [];
    allSims = [];
    trainSet = bsxfun(@plus,[0:(setSize-1)],[1:stepSize:length(sim)-(setSize-1)]');
    for ti = 1:length(trainSet(:,1))
        test = setxor([1:length(sim)],trainSet(ti,:));
        [simVals guess] = nanmax(sim(trainSet(ti,:),test),[],1);
        trainEnvs = envs(trainSet(ti,:));
        predictions = trainEnvs(guess);
        actual = envs(test);
        isEnv = strcmp(predictions,actual)';
        isContextGroup = strcmp(cellfun(@(x) x(1),predictions,'un',0), ...
            cellfun(@(x) x(1),actual,'un',0))';
        
        lags = nanmin(abs(bsxfun(@minus,trainSet(ti,:)',test)),[],1);
        allEnvAcc = [allEnvAcc [lags; isEnv]];
        allCGAcc = [allCGAcc [lags; isContextGroup]];
        allSims = [allSims [lags; simVals]];
    end
    
    blah = [1:6:nanmax(allEnvAcc(1,:))];
    blah(end) = inf;
    lagGroup = [blah(1:end-1)' blah(2:end)'-1];
    accuracies = nan(length(lagGroup(:,1)),3);
    for i = 1:length(lagGroup(:,1))
        accuracies(i,1) = nanmean(allEnvAcc(2,allEnvAcc(1,:)>=lagGroup(i,1) & ...
            allEnvAcc(1,:)<=lagGroup(i,2)));
        accuracies(i,2) = nanmean(allCGAcc(2,allCGAcc(1,:)>=lagGroup(i,1) & ...
            allCGAcc(1,:)<=lagGroup(i,2)));
        accuracies(i,3) = nanmean(allSims(2,allSims(1,:)>=lagGroup(i,1) & ...
            allSims(1,:)<=lagGroup(i,2)));
    end
end