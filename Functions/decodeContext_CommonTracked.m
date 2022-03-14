function accuracies = decodeContext_CommonTracked(aSamp,uP,uGT,envs,stepSize,setSize)
    tic
    if nargin < 3 || isempty(stepSize)
        stepSize = 6;
    end
    if nargin < 4 || isempty(setSize)
        setSize = 6;
    end
    strLength = 0;
    fprintf('\n\tDecoding from commonly-tracked cells (matched sampling)...  ')
    
    allEnvAcc = [];
    allCGAcc = [];
    trainSet = bsxfun(@plus,[0:(setSize-1)],[1:stepSize:32-(setSize-1)]');
    for ti = 1:length(trainSet(:,1))
        
        fprintf(repmat('\b',[1 strLength]));
        str = sprintf(['(Training Set:  ' num2str(ti) ')']);
        fprintf(str);
        strLength = length(str);
        
        test = setxor([1:32],trainSet(ti,:));
        
        tmp = cellfun(@all,cellfun(@isnan,uGT,'un',0),repmat({2},size(uGT)),'un',0);
        activeCells = ~cat(2,tmp{:});
        includeCells = bsxfun(@and,all(activeCells(:,trainSet(ti,:)),2),activeCells(:,test));
        trainEnvs = envs(trainSet(ti,:));
        predictions = [];
        for i = 1:length(test)
            vals = nan(1,setSize);
            for k = 1:length(trainSet(1,:))
                
                doSamps = round(nanmin(aSamp(:,:,[trainSet(ti,k) test(i)]),[],3).*30);
                testMap = mkTraceMaps(uP{test(i)}, ...
                        uGT{test(i)}(includeCells(:,i),:),[],[17 17],doSamps,1);
                rt1 = reshape(testMap,[numel(testMap(:,:,1)) length(testMap(1,1,:))]);
                
                trainMap = mkTraceMaps(uP{trainSet(ti,k)}, ...
                    uGT{trainSet(ti,k)}(includeCells(:,i),:),[],[17 17],doSamps,1);
                rt2 = reshape(trainMap,[numel(trainMap(:,:,1,1)) length(trainMap(1,1,:,1))]);
                
                isGoodPixels = ~all(isnan(rt1),2) & ~all(isnan(rt2),2);
                
                art1 = rt1(isGoodPixels,:);
                art2 = rt2(isGoodPixels,:);
                vals(k) = corr(art1(:),art2(:));
            end
            [a b] = nanmax(vals);
            predictions = [predictions; trainEnvs(b)];
        end
        actual = envs(test);
        isEnv = strcmp(predictions,actual)';
        isContextGroup = strcmp(cellfun(@(x) x(1),predictions,'un',0), ...
            cellfun(@(x) x(1),actual,'un',0))';
        
        lags = nanmin(abs(bsxfun(@minus,trainSet(ti,:)',test)),[],1);
        allEnvAcc = [allEnvAcc [lags; isEnv]];
        allCGAcc = [allCGAcc [lags; isContextGroup]];
    end
    
    lagGroup = [1 6; 7 12; 13 18; 19 inf];
    accuracies = nan(length(lagGroup(:,1)),2);
    for i = 1:length(lagGroup(:,1))
        accuracies(i,1) = nanmean(allEnvAcc(2,allEnvAcc(1,:)>=lagGroup(i,1) & ...
            allEnvAcc(1,:)<=lagGroup(i,2)));
        accuracies(i,2) = nanmean(allCGAcc(2,allCGAcc(1,:)>=lagGroup(i,1) & ...
            allCGAcc(1,:)<=lagGroup(i,2)));
    end
    
    tmp = toc;
    fprintf('  Time:  %0.3fs.',tmp);
end