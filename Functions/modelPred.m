function modelPred()
    clc
    close all
    drawnow
    
%     train = [1:12]; % 27:36 [1:14]  [41:54]
    for doTrain = fliplr([{[1:12]} {[42:53]} {[96-11:96]}]);
        train = doTrain{1};
        cGroup = [{[{'sq1'} {'sq2'} {'sq3'}]} {[{'g1'} {'g2'} {'g3'}]}];

        
        
        paths = getFilePaths('MatlabData/Model','.mat');
        simGroup = [];
    
    
        doPixels = 2:3:48;
        predAcc = [];
        fprintf('\t\t\tPredicting context from simulated data:\n')
        for p = paths'
            fprintf(['\n\t' p{1}])
            load(p{1},'ca1Maps','envs','envLabels');

            fn = p{1}(find(ismember(p{1},'/'),1,'last')+1:end-4);
            if all(fn(1:7)=='shuffle')
                simGroup = [simGroup; 1];
            else
                simGroup = [simGroup; 2];
            end

            ndays = 1:length(ca1Maps(1,1,1,:));

            % Setup training set

            tr = ca1Maps(doPixels,doPixels,:,train).*100;
            tr = reshape(tr,[numel(tr(:,:,1,1)) length(tr(1,1,:,1)) length(train)]);
    %         tr = poissrnd(tr);
            tmp = permute(tr,[1 3 2]);
            trVecs = reshape(tmp,[numel(tmp(:,:,1)) length(tmp(1,1,:))]);

            isGroup = zeros(1,length(train));
            trenvs = envs(train);
            for i = 1:length(trenvs)
                if ismember(trenvs{i},cGroup{1})
                    isGroup(i) = 1;
                end
                if ismember(trenvs{i},cGroup{2})
                    isGroup(i) = 2;
                end
            end

            blah = repmat(isGroup,[numel(ca1Maps(doPixels,doPixels,1,1)) 1]);
            traingroup = blah(:);

            % set up testing set

            test = ndays(~ismember(ndays,train));

            tr = ca1Maps(doPixels,doPixels,:,test).*100;
            tr = reshape(tr,[numel(tr(:,:,1,1)) length(tr(1,1,:,1)) length(test)]);
    %         tr = poissrnd(tr);
            tmp = permute(tr,[1 3 2]);
            testVecs = reshape(tmp,[numel(tmp(:,:,1)) length(tmp(1,1,:))]);

            isGroup = zeros(1,length(test));
            trenvs = envs(test);
            for i = 1:length(trenvs)
                if ismember(trenvs{i},cGroup{1})
                    isGroup(i) = 1;
                end
                if ismember(trenvs{i},cGroup{2})
                    isGroup(i) = 2;
                end
            end

            blah = repmat(isGroup,[numel(ca1Maps(doPixels,doPixels,1,1)) 1]);
            testgroup = blah(:);

            %%%%%% Do svm things
            
            %eliminate nan pixels from training
            exclude = all(isnan(trVecs),2);
            trVecs(exclude,:) = [];
            traingroup(exclude) = [];
            
            nsims = 100;
            allPredXtime = nan(nsims,length(test));
            for q = 1:nsims

                trainNum = 50;

                doTrain1 = find(traingroup==1);
                doTrain2 = find(traingroup==2);

                t1 = randperm(length(doTrain1));
                t2 = randperm(length(doTrain2));

                doTrain = [doTrain1(t1(1:trainNum)); doTrain2(t2(1:trainNum))];

                outSVM = svmtrain(trVecs(doTrain,:),traingroup(doTrain,:));

                predGroup = svmclassify(outSVM,testVecs);

                correctPrediction = predGroup == testgroup;
                exclude = reshape(all(isnan(testVecs),2),[numel(ca1Maps(doPixels,doPixels,1,1)) length(test)]);
                tmp = double(reshape(correctPrediction,[numel(ca1Maps(doPixels,doPixels,1,1)) length(test)]));
                tmp(exclude) = nan;
                predXtime = nanmean(tmp);
                allPredXtime(q,:) = predXtime;
            end

            figure(1)
            set(gcf,'position',[50 50 400 250])
            plot(test,nanmean(allPredXtime));
            set(gca,'ylim',[0 1]);
            hold on
            plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--')
            drawnow


            predAcc = [predAcc; nanmean(allPredXtime)];
        end

        close all

        low = predAcc(simGroup==1,:);
        high = predAcc(simGroup==2,:);
        figure(1)
        set(gcf,'position',[50 50 400 250])
        h = mkLine([{low(1:80,:)} {high(1:80,:)}],test);
        set(gca,'ylim',[0 1],'xlim',[nanmin(ndays)-1 nanmax(ndays)+1],'xtick',[8:8:96],'xticklabel',[8:8:96]);
        hold on
        h2 = patch([train(1) train(end) train(end) train(1)],[0 0 1 1],[0.5 0.5 0.5],...
            'facecolor',[0.5 0.5 0.5],'edgecolor','none');
        plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--')
        legend(h,[{'Low selectivity'} {'High selectivity'}],'location','southeast')    
        pvals = nan(1,length(predAcc(1,:)));
        for k = 1:length(predAcc(1,:))
            [a b c] = ranksum(low(1:80,k),high(1:80,k));
            pvals(k) = a;
        end
        text(test(pvals<0.05),0.9.*ones(1,nansum(pvals<0.05)),'*','fontname','arial', ...
            'fontweight','normal','fontsize',9,'horizontalalignment','center')
        xlabel('Day')
        ylabel('Prediction accuracy (%)')
        set(gca,'ytick',[0:0.25:1],'yticklabels',[0:25:100])
        drawnow

        outP = ['Plots/Model/Summary/Decoding_Train_' num2str(train(1)) '_' num2str(train(2))];
        saveFig(gcf,outP,[{'tiff'} {'pdf'}])
    end
end





























