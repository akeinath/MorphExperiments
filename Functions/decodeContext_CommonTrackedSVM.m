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
        
        isGroup = strcmp(cellfun(@(x) x(1),envs,'un',0),{'g'})';
        
        test = setxor([1:32],trainSet(ti,:));
        
        tmp = cellfun(@all,cellfun(@isnan,uGT,'un',0),repmat({2},size(uGT)),'un',0);
        activeCells = ~cat(2,tmp{:});
        includeCells = bsxfun(@and,all(activeCells(:,trainSet(ti,:)),2),activeCells(:,test));
        trainEnvs = envs(trainSet(ti,:));
        accuracy = nan(1,length(test));
        for i = 1:length(test)
            vals = nan(1,setSize);
            doSamps = round(nanmin(aSamp(:,:,[trainSet(ti,:) test(i)]),[],3).*30);
            testMap = mkTraceMaps(uP{test(i)}, ...
                    uGT{test(i)}(includeCells(:,i),:),[],[17 17],doSamps,1);
            rt1 = reshape(testMap,[numel(testMap(:,:,1)) length(testMap(1,1,:))]);
            
            trainMaps = nan([size(testMap) length(trainSet(1,:))]);
            for k = 1:length(trainSet(1,:))                
                trainMaps(:,:,:,k) = mkTraceMaps(uP{trainSet(ti,k)}, ...
                    uGT{trainSet(ti,k)}(includeCells(:,i),:),[],[17 17],doSamps,1);
            end
            
            rt2 = reshape(trainMaps,[numel(trainMaps(:,:,1,1)) ...
                length(trainMaps(1,1,:,1)) length(trainMaps(1,1,1,:))]);

            isGoodPixels = ~all(isnan(rt1),2);

            art1 = rt1(isGoodPixels,:);
            art2 = rt2(isGoodPixels,:,:);
            vals = nan(1,length(trainSet(1,:)));
            for k = 1:length(trainSet(1,:))
                tmp = art2(:,:,k);
                vals(k) = corr(art1(:),tmp(:));
            end
            
            tart2 = reshape(permute(art2,[1 3 2]),[length(art2(:,1,1)).*length(art2(1,1,:)) ...
                length(art2(1,:,1))]);

            tmp = repmat(isGroup(trainSet(ti,:)),[nansum(isGoodPixels) 1]);
            trLabels = tmp(:);
            
            outSVM = fitcsvm(tart2,trLabels(:));
            pred = predict(outSVM,art1);
            accuracy(i) = nanmean(pred==isGroup(test(i)));
        end
        lags = nanmin(abs(bsxfun(@minus,trainSet(ti,:)',test)),[],1);
        allCGAcc = [allCGAcc [lags; accuracy]];
    end
    
    blah = [1:1:nanmax(allCGAcc(1,:))];
    blah(end) = inf;
    lagGroup = [blah(1:end-1)' blah(2:end)'-1];
    accuracies = nan(length(lagGroup(:,1)),3);
    for i = 1:length(lagGroup(:,1))
        accuracies(i,2) = nanmean(allCGAcc(2,allCGAcc(1,:)>=lagGroup(i,1) & ...
            allCGAcc(1,:)<=lagGroup(i,2)));
    end
    
    tmp = toc;
    fprintf('  Time:  %0.3fs.',tmp);
end