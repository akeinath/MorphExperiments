function [correct correctSim] = help_matchAndPredict(train_vec,train_group,test_vec,test_group,nsims,useCellNum)
        
    if nargin<5 || isempty(nsims)
        nsims = 1;
    end
    
    if nargin<6 || isempty(useCellNum)
        useCellNum = 1;
    end
    
    correct = nan([nsims 1]);
    correctSim = nan([nsims 1]);
%     meanAbsoluteWeights = 1;
    for q = 1:nsims

        doTrain = repmat({[]},[1 nanmax(train_group(:))]);
        t = repmat({[]},[1 nanmax(train_group(:))]);
        for k = 1:nanmax(train_group(:))
            doTrain{k} = find(train_group==k);
            t{k} = randperm(length(doTrain{k}));
        end

        trainNum = nanmin(cellfun(@length,doTrain));

        tmp = [];
        for k = 1:nanmax(train_group(:))
            tmp = [tmp; doTrain{k}(t{k}(1:trainNum))];
        end
        doTrain = tmp;

% % %         hiddenLayerSize = 20;
% % %         net = patternnet(hiddenLayerSize);
% % %         % Set up Division of Data for Training, Validation, Testing
% % %         net.divideParam.trainRatio = 0.5;
% % %         net.divideParam.valRatio = 0.5;
% % %         net.divideParam.testRatio = 0;
% % %         net.trainParam.showWindow = 0;
% % %         net.biasConnect = false(2,1);
% % %         [net,tr] = train(net,train_vec(:,doTrain(:)),train_group(doTrain(:))-1);
% % %         outputs = net(train_vec(:,doTrain));
% % % %         performance = perform(net,train_group(doTrain(:))'-1,outputs);
% % %         outputs = [net(test_vec)>0.5]+1;
% % %         correct(q) = (outputs == test_group');
% % % %             nanmean(abs(net.IW{1}));

        if ~isempty(train_vec) && length(train_vec(:,1))>1
            xc = corr(train_vec(:,doTrain(:)),test_vec);
            tmpTrainGroup = train_group(doTrain(:));
            [a b] = nanmax(xc,[],1);
            
            tac = nan(1,length(test_group));
            for i = 1:length(test_group)
                tac(i) = tmpTrainGroup(b(i)) == test_group(i);
            end
            correct(q) = nanmean(tac);
%             correct(q) = tmpTrainGroup(b) == test_group;
            correctSim = a;
        end
    end
%     mkBow([{nnc} {correct}])
    correct = nanmean(correct);
    correctSim = nanmean(correctSim);
end


% % %                 doAgain = true;
% % %                 while doAgain
% % %                 
% % %                     hiddenLayerSize = 1;
% % %                     net = patternnet(hiddenLayerSize);
% % %                     % Set up Division of Data for Training, Validation, Testing
% % %                     net.divideParam.trainRatio = 0.8;
% % %                     net.divideParam.valRatio = 0.2;
% % %                     net.divideParam.testRatio = 0;
% % %                     net.trainParam.showWindow = 0;
% % %                     net.biasConnect = false(2,1);
% % %                     [net,tr] = train(net,train_vec(doTrain,doCells)',train_group(doTrain)'-1);
% % %                     outputs = net(train_vec(doTrain,doCells)');
% % %                     performance = perform(net,train_group(doTrain)'-1,outputs);
% % %                     if performance < 0.01
% % %                         doAgain = false;
% % %                     end
% % %                 end
% % %                 outputs = [net(test_vec(:,doCells)')>0.5]+1;
% % %                 accXday(q,testI) = nanmean(outputs==test_group');
                
                
%                 outModel = fitcsvm(train_vec(doTrain,doCells),train_group(doTrain)-1, ...
%                     'Standardize',false,'KernelFunction','RBF','KernelScale','auto');
%                 predGroup = predict(outModel,test_vec(:,doCells));
%                 accXday(q,testI) = nanmean(predGroup == (test_group-1));
                
%                 xc = pdist2(train_vec(doTrain,doCells),test_vec(:,doCells),'cosine');