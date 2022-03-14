function modelPred()
    clc
    close all
    drawnow
    
%     train = [1:12]; % 27:36 [1:14]  [41:54]
    cGroup = [{[{'sq1'} {'sq2'} {'sq3'}]} {[{'g1'} {'g2'} {'g3'}]}];

    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve

    paths = getFilePaths('MatlabData/Model_Grouped_SVMs_12Days','.mat');
%     paths(1:80) = [];
    
    doPixels = 2:3:48;
    allPredAcc = [];
    meanSel = [];
    propHighSel = [];
    fprintf('\t\t\tPredicting context from simulated data:\n')
    for p = paths'
        close all
        drawnow
        fprintf(['\n\t' p{1}])
        
        outP = ['MatlabData/Model_Grouped_SVMs_18Days/' p{1}(find(ismember(p{1},'/'),1,'last')+1:end-4)];
%         if exist([outP '.mat'])==2
%             fprintf('\n\t\tAlready Ran');
%             continue
%         end
        
        load(p{1},'ca1Maps','envs','envLabels','weights');


        %%%%%%%%%%%%%%%%%% Get Cell Measures %%%%%%%%%%%%%%%%%%%%%%%

%          [r2penalty tr2] = help_Maps2Var(ca1Maps,envs,envLabels);

        inds = slind(p{1});
        root = ['Plots/' p{1}(inds(1)+1:end-4)];
        
        useSessions = 1:32;
        sim = getPairwiseMapSim(ca1Maps(:,:,:,useSessions));
        mds2D(sim.pearson,envs(useSessions),envLabels,[root '/MDS_2D']);
        
%          tmp = sort(r2penalty,2);
%          selectivity = [tmp(:,3)-nanmean(tmp(:,1:2),2)];
%          meanSel = [meanSel; nanmean(selectivity)];
%          propHighSel = [propHighSel; nanmean(selectivity>0.4)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %         ndays = 1:length(ca1Maps(1,1,1,:));
% % % 
% % %         % Setup training set
% % % 
% % %         allVecs = ca1Maps(doPixels,doPixels,:,:);
% % %         allVecs = reshape(allVecs,[numel(allVecs(:,:,1,1)) length(allVecs(1,1,:,1)) length(ndays)]);
% % % 
% % %         isGroup = zeros(1,length(ndays));
% % %         trenvs = envs(:);
% % %         for i = 1:length(trenvs)
% % %             if ismember(trenvs{i},cGroup{1})
% % %                 isGroup(i) = 1;
% % %             end
% % %             if ismember(trenvs{i},cGroup{2})
% % %                 isGroup(i) = 2;
% % %             end
% % %         end
% % % 
% % %         [predAcc linAcc] = help_slidingSVM(allVecs,isGroup);
% % %         
% % %         checkP(outP);
% % %         save(outP);
% % %         
% % % %         f = @(b,x) b(1).*exp(b(2).*x)+b(3);                                     % Objective Function
% % % %         B = fminsearch(@(b) nanmean((apred - f(b,test)).^2), [1; -1; 0])     
% % % 
% % % %             [params fval] = fmincon(@(p)exp_error(test,apred,p),[0 1 1] ,[],[],[],[],...
% % % %                 [-inf 1e-10 1e-10],[inf inf inf],[],optimoptions('fmincon','Display','none'));
% % % % 
% % % %         figure(1)
% % % %         set(gcf,'position',[50 50 400 250])
% % % %         plot(test,apred,'color','k');
% % % %         set(gca,'ylim',[0 1]);
% % % %         hold on
% % % %         plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--')
% % % %         plot(test,f(B,test),'color','r','linestyle','-')
% % % %         drawnow
    end
end

function error = exp_error(test,data,p)
    pred = doExp(test,p);
    error = nanmean((data-pred).^2);
end

function pred = doExp(x,p)
    pred =  -p(1)+p(2).*exp(-p(3)*x);
end




























