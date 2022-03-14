function modelPred()
    clc
    close all
    drawnow
    
%     train = [1:12]; % 27:36 [1:14]  [41:54]
    cGroup = [{[{'sq1'} {'sq2'} {'sq3'}]} {[{'g1'} {'g2'} {'g3'}]}];

    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve

    for blah = [{'UniformDrift'} {'EnvSpecificDrift'} {'AttractGroupDrift'}]
    
        paths = getFilePaths(['MatlabData/Modeling_GroupedDrift/' blah{1}],'.mat');
    %     paths(1:80) = [];

        doPixels = 2:3:48;
        allPredAcc = [];
        allPredSim = [];
        fprintf('\t\t\tPredicting context from simulated data:\n')
        allAngles = [];
        allDriftVariance = [];
        for p = paths'
            close all
            drawnow
            fprintf(['\n\t' p{1}])

            outP = ['MatlabData/Model_DifferentDrift/' p{1}(find(ismember(p{1},'/'),1,'last')+1:end-4)];
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
            [a b] = mds2D(sim.pearson,envs(useSessions),envLabels,[root '/MDS_2D'],true);
            allAngles = [allAngles; a];
            allDriftVariance = [allDriftVariance; b];


    % % %         [accXset accXset_sim] = help_contextClassification(ca1Maps,envs,[33:38],false,false);
    % % %         allPredAcc = [allPredAcc; accXset];
    % % %         allPredSim = [allPredSim; accXset_sim];

    %         close all
    %         figure
    %         set(gcf,'position',[50 50 700 250])
    %         subplot(1,2,1)
    %         mkLine(allPredAcc);
    %         subplot(1,2,2)
    %         mkLine(allPredSim);
    %         drawnow
        end

        save([blah{1} '_Angles'],'allAngles','allDriftVariance');
        
    end
    
    
    
    a = load('EnvSpecificDrift_Angles');
    b = load('AttractGroupDrift_Angles');
    c = load('UniformDrift_Angles');
    
%     mkPolar([{mod(a.allDriftVariance(1:50),180).*2} ...
%         {mod(b.allDriftVariance(1:50),180).*2} ...
%         {mod(c.allDriftVariance(1:50),180).*2}],10)

    tmpA = nanmin(abs(a.allDriftVariance),360-abs(a.allDriftVariance));
    tmpB = nanmin(abs(b.allDriftVariance),360-abs(b.allDriftVariance));
    tmpC = nanmin(abs(c.allDriftVariance),360-abs(c.allDriftVariance));
    
%     figure(1)
%     [h at] = groupHist([{tmpA} {tmpB} {tmpC}],[0:5:180])
    
    mkPolar([{a.allDriftVariance(1:50)} ...
        {b.allDriftVariance(1:50)} ...
        {c.allDriftVariance(1:50)}],10)
    inds = slind(p{1});
    root = ['Plots/' p{1}(inds(1)+1:inds(2)-1)];
    saveFig(gcf,[root '/Summary/MDS_AngleDistributions'],[{'tiff'} {'pdf'} {'jpeg'}]);
end

function error = exp_error(test,data,p)
    pred = doExp(test,p);
    error = nanmean((data-pred).^2);
end

function pred = doExp(x,p)
    pred =  -p(1)+p(2).*exp(-p(3)*x);
end




























