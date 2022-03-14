function [accXset accXset_sim] = help_contextClassification(maps,envs,trainSets,shuffleTest,crossValSelection,cGroup)

    if nargin < 4 || isempty(shuffleTest)
        shuffleTest = false;
    end
    
    if nargin < 5 || isempty(crossValSelection)
        crossValSelection = false;
    end
    
    if nargin < 6 || isempty(cGroup)
        cGroup = [{[{'sq1'} {'sq2'} {'sq3'}]} {[{'g1'} {'g2'} {'g3'}]}];
    end
    
    fprintf('\n\t\tPredicting context group from data:\t')
    accXset = nan([length(trainSets(:,1)) length(maps(1,1,1,:))-length(trainSets(1,:))]);
    accXset_sim = nan([length(trainSets(:,1)) length(maps(1,1,1,:))-length(trainSets(1,:))]);
    clearLen = 0;
    
    activeCells = (~permute(all(all(isnan(maps),1),2),[3 4 1 2]));

    isGroup = zeros(length(envs),1);
    for i = 1:length(cGroup)
        isGroup(ismember(envs,cGroup{i})) = i;
    end
    
    for trainVecI = 1:length(trainSets(:,1))
        str = sprintf([num2str(trainVecI) ' of ' num2str(length(trainSets(:,1)))]);
        fprintf([repmat('\b',[1 clearLen]) str]);
        clearLen = length(str);
        
        train1 = trainSets(trainVecI,:);
                
        ndays = 1:length(maps(1,1,1,:));

        
%         [blah isGroup] = ismember(envs,unique(envs)');
        
        [train_vec train_group train_cellInds] = help_mapSVMPrep(maps,isGroup,train1);

        if crossValSelection
            [cross_vec_1 cross_1_group cross_1_cellInds] = help_mapSVMPrep(maps,isGroup,train1(1:round(length(train1)./2)));
            [cross_vec_2 cross_2_group cross_2_cellInds] = help_mapSVMPrep(maps,isGroup,train1(round(length(train1)./2)+1:end));
        
            stillActCells = find(all(activeCells(:,train1),2));
            
            allCellVals = nan(length(stillActCells),length(stillActCells));
            for numCells = 10 % 5:5:length(stillActCells)
                crossSims = 10000;
                cross_ensembles = nan(crossSims,numCells);
                cross_acc = nan(crossSims,1);
                for i = 1:crossSims
                    tmp = stillActCells(randperm(length(stillActCells)));
                    
                    ttrv = cross_vec_1(ismember(cross_1_cellInds,tmp(1:numCells)'),:);
                    ttv = cross_vec_2(ismember(cross_2_cellInds,tmp(1:numCells)'),:);

                    excludePixels = any(isnan(ttrv),2)|any(isnan(ttv),2);
                    ttrv(excludePixels,:) = [];
                    ttv(excludePixels,:) = [];
                    
                    cross_acc(i) = help_matchAndPredict(ttrv, ...
                        cross_1_group,ttv,cross_2_group,1);
                    cross_ensembles(i,:) = tmp(1:numCells);
                end
                
                
                cellVals = nan(length(stillActCells),1);
                for k = 1:length(stillActCells)
                    cellVals(k,:) = nanmean(cross_acc(any(cross_ensembles==stillActCells(k),2)));
                end
                
                allCellVals(:,numCells) = cellVals;
            end
            cellValRank = normRank(allCellVals(:,10))./length(allCellVals(:,1));
            crossValActive = stillActCells;
        end
        
        % set up testing set
        inc = 0.5; % crossval increment
        test = ndays(~ismember(ndays,train1));
        correct = nan(1,length(test));
        correct_sim = nan(1,length(test));
        
        if crossValSelection
            correct = nan(round(1./inc),length(test));
            correct_sim = nan(round(1./inc),length(test));
        end
        
        if shuffleTest
            tmp = isGroup(test);
            isGroup(test) = tmp(randperm(length(tmp)));
        end
        
        for testI = 1:length(test)
            
            [test_vec test_group test_cellInds] = help_mapSVMPrep(maps,isGroup,test(testI));   
            
%             if ~crossValSelection
                stillActCells = find(all(activeCells(:,[train1 test(testI)]),2));
                ttrv = train_vec(ismember(train_cellInds,stillActCells'),:);
                ttv = test_vec(ismember(test_cellInds,stillActCells'),:);

                excludePixels = any(isnan(ttrv),2)|any(isnan(ttv),2);
                ttrv(excludePixels,:) = [];
                ttv(excludePixels,:) = [];

                [a b] = help_matchAndPredict(ttrv, ...
                    train_group,ttv,test_group,1);
                correct(testI) = a;
                correct_sim(testI) = b;
%             else
%                 stillActCells = find(all(activeCells(:,[train1 test(testI)]),2));
%                 for rv = inc:inc:1
%                 
%                     doCells = stillActCells(ismember(stillActCells, ...
%                         crossValActive(cellValRank>=rv-inc & cellValRank<rv)));
%                     
%                     ttrv = train_vec(ismember(train_cellInds,doCells'),:);
%                     ttv = test_vec(ismember(test_cellInds,doCells'),:);
% 
%                     excludePixels = any(isnan(ttrv),2)|any(isnan(ttv),2);
%                     ttrv(excludePixels,:) = [];
%                     ttv(excludePixels,:) = [];
% 
%                     [a b] = help_matchAndPredict(ttrv, ...
%                         train_group,ttv,test_group,1);
%                     
%                     correct(round(rv./inc),testI) = a;
%                     correct_sim(round(rv./inc),testI) = b;
%                 end
%             end
        end
        accXset(trainVecI,:) = correct;
        accXset_sim(trainVecI,:) = correct_sim;
    end
    
    tmp = [nan(size(accXset)) accXset];
    for i = 1:length(tmp(:,1))
        tmp(i,:) = circshift(tmp(i,:),[0 -(i-1)]);
    end
    tmp = cat(3,tmp(:,length(accXset(1,:))+1:end),fliplr(tmp(:,1:length(accXset(1,:)))));
    accXset = (nanmean(reshape(permute(tmp,[1 3 2]),[numel(tmp(:,1,:)) length(tmp(1,:,1))]),1));
    
    
    tmp = [nan(size(accXset_sim)) accXset_sim];
    for i = 1:length(tmp(:,1))
        tmp(i,:) = circshift(tmp(i,:),[0 -(i-1)]);
    end
    tmp = cat(3,tmp(:,length(accXset_sim(1,:))+1:end),fliplr(tmp(:,1:length(accXset_sim(1,:)))));
    accXset_sim = (nanmean(reshape(permute(tmp,[1 3 2]),[numel(tmp(:,1,:)) length(tmp(1,:,1))]),1));
end

































































