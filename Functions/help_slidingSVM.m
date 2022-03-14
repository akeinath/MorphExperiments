function [predAcc linAcc] = help_slidingSVM(allVecs,group)

    group = repmat(group,[length(allVecs(:,1,1)) 1 1]);
    day = repmat(1:length(allVecs(1,1,:)),[length(allVecs(:,1,1)) 1 1]);
    group = group(:);
    day = day(:);
    
    allVecs= permute(allVecs,[1 3 2]);
    allVecs = reshape(allVecs,[numel(allVecs(:,:,1)) length(allVecs(1,1,:))]);
    
    exclude = all(isnan(allVecs),2);
    allVecs(exclude,:) = [];
    group(exclude) = [];
    day(exclude) = [];
    
    tic
    windowSize = 12;
    predAcc = nan(nanmax(day)-windowSize+1,nanmax(day));
    for train = 1:nanmax(day)-windowSize+1
        
        isTrain = day>=train & day <train+windowSize;
        
        doTrain1 = find(group==1 & isTrain);
        doTrain2 = find(group==2 & isTrain);
        
        nsims = 30;
        allPredXtime = nan(nsims,nanmax(day));
        for q = 1:nsims

            trainNum = 50;

            t1 = randperm(length(doTrain1));
            t2 = randperm(length(doTrain2));

            doTrain = [doTrain1(t1(1:trainNum)); doTrain2(t2(1:trainNum))];

            outSVM = svmtrain(allVecs(doTrain,:),group(doTrain,:));

            predGroup = svmclassify(outSVM,allVecs(~isTrain,:));

            correctPrediction = predGroup == group(~isTrain);
            td = day(~isTrain);
            acc = nan(1,nanmax(day));
            for k = 1:nanmax(day)
                acc(k) = nanmean(correctPrediction(td==k));
            end
            
            allPredXtime(q,:) = acc;
        end
        
        apred = nanmean(allPredXtime,1);
        predAcc(train,:) = apred;
    end
    cspredAcc = [nan(size(predAcc)) predAcc];
    for i = 1:length(cspredAcc(:,1))
        cspredAcc(i,:) = circshift(cspredAcc(i,:),[0 -(i-1)]);
    end
    toc
    
    a = cspredAcc(:,length(predAcc(1,:))+1:end);
    b = (circshift(fliplr(cspredAcc(:,1:length(predAcc(1,:)))),[0 windowSize]));
    linAcc = circshift(nanmean([a; b]),[0 -windowSize]);
end