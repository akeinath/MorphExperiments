function accuracies = decodeContext_MapSVM(um,envs,stepSize,setSize)
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
    trainSet = bsxfun(@plus,[0:(setSize-1)],[1:stepSize:length(um(1,1,1,:))-(setSize-1)]');
    for ti = 1:length(trainSet(:,1))
        
        fprintf(repmat('\b',[1 strLength]));
        str = sprintf(['(Training Set:  ' num2str(ti) ')']);
        fprintf(str);
        strLength = length(str);
        
        isGroup = strcmp(cellfun(@(x) x(1),envs,'un',0),{'g'})';
        
        test = setxor([1:length(um(1,1,1,:))],trainSet(ti,:));
        
        trainMaps = um(:,:,:,trainSet(ti,:));
        rt2 = reshape(trainMaps,[numel(trainMaps(:,:,1,1)) ...
            length(trainMaps(1,1,:,1)) length(trainMaps(1,1,1,:))]);
        isGoodPixels = all(~all(isnan(rt2),2),3);
        rt1 = reshape(um(:,:,:,test),[numel(um(:,:,1,1)) length(um(1,1,:,1)) length(test)]);
        
        art1 = rt1(isGoodPixels,:,:);
        art2 = rt2(isGoodPixels,:,:);
        
        tart2 = reshape(permute(art2,[1 3 2]),[length(art2(:,1,1)).*length(art2(1,1,:)) ...
                length(art2(1,:,1))]);
        tmp = repmat(isGroup(trainSet(ti,:)),[nansum(isGoodPixels) 1]);
        trLabels = tmp(:);
        outSVM = fitcsvm(tart2,trLabels(:));
        
        tart1 = reshape(permute(art1,[1 3 2]),[length(art1(:,1,1)).*length(art1(1,1,:)) ...
                length(art1(1,:,1))]);
        tmp = repmat(isGroup(test),[nansum(isGoodPixels) 1]);
        testLabels = tmp(:);

        pred = predict(outSVM,tart1);
        accuracy = nanmean(reshape(pred==testLabels,[length(art1(:,1,1)) length(test)]),1);
        
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