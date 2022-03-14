function batchAnalyses_Upload(paths)

    if nargin < 1 || isempty(paths)
        paths = getFilePaths('MatlabData/AnalysisPreloads/MorphParadigm_CA1_RiseExtracted','.mat');
    end

    across = prepAcross();
    
    warning off all
    envLabel = [{'sq1'} {'sq2'} {'sq3'} {'g3'} {'g2'} {'g1'}];
    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve
    
    clc
    fprintf('\n\t\t\t*********Running batched analyses********\n\n');
    for mi = 1:length(paths)
        close all
        drawnow
        fprintf(['\n' paths{mi}]);
        load(paths{mi},'um','umfr','upfr','ushc','envs','uIsPC','SFPs','aSamp', ...
            'matchedSim','matchedBothSim','aSIC','uP'); % 'uGT','uP'
        slashInds = find(ismember(paths{mi},'/'));
        root = ['Plots/BatchedAnalyses/' paths{mi}(slashInds(3)+1:end-4)];
        
        across.maps = cat(3,across.maps,um);
        across.envs = cat(1,across.envs,repmat(envs',[size(um,3) 1]));
        across.mfr = cat(1,across.mfr,umfr);
                
        % Sampling stuff
%         batchPlotPaths(uP,[root '/Paths']);
        across.coverage = [across.coverage; relativeCoverage(aSamp,1,envs)];
        across.isPlaceCell = [across.isPlaceCell; uIsPC];

        % MDS stuff

        [angdiff driftVariance driftAmount stress] = mds2D( ...
            matchedBothSim.pearson,envs,envLabel,[root '/MDS_2D_matchedSampling_and_Comparisons']);
        across.mds.angles = [across.mds.angles; angdiff];
        across.mds.driftVariability = [across.mds.driftVariability; driftVariance];
        across.mds.driftAmount = cat(1,across.mds.driftAmount,driftAmount'); 
        across.mds.stress = cat(1,across.mds.stress,stress); 
        mds3D(matchedBothSim.pearson,envs,envLabel,[root '/MDS_3D_matchedSampling_and_Comparisons']);
        [angdiff driftVariance driftAmount stress] = mds2D( ...
            matchedBothSim.pv,envs,envLabel,[root '/MDS_PV/MDS_2D_matchedSampling_and_Comparisons']);
        mds3D(matchedBothSim.pv,envs,envLabel,[root '/MDS_PV/MDS_3D_matchedSampling_and_Comparisons']);
        
        [transitionPoint maxDecorr] = sim2seq(matchedSim.pearson,uIsPC<0.05,envs,doComps,root);
        across.seqAnalysis.maxDecorr = [across.seqAnalysis.maxDecorr; maxDecorr];
        across.seqAnalysis.transitionPoint = [across.seqAnalysis.transitionPoint; transitionPoint];

        % GLM Analysis
        % 2D
        minSessions = 16;
        morphPoint = ceil(transitionPoint);
        morphPoint(end+1) = morphPoint(end);
        RDMs = getRDMs_r2(envs,envLabel,morphPoint);
        gRDMs = [{-RDMs(:,:,1:3)} {-RDMs(:,:,4:end)}];
        [r2penalty tr2 shuffle_r2penalty shuffle_tr2] = classFactorLoading( ...
            matchedSim.pearson,gRDMs,nchoosek(minSessions,2));
        across.class.d2.r2penalty = [across.class.d2.r2penalty; r2penalty];
        across.class.d2.tr2 = [across.class.d2.tr2; tr2];
        across.class.d2.shuffle_r2penalty = [across.class.d2.shuffle_r2penalty; shuffle_r2penalty];
        across.class.d2.shuffle_tr2 = [across.class.d2.shuffle_tr2; shuffle_tr2];
        plotgRDMs(gRDMs,[root '/RDMs_2D']);

        % 3D
        minSessions = 16;
        morphPoint = ceil(transitionPoint);
        morphPoint(end+1) = morphPoint(end);
        RDMs = getRDMs(envs,envLabel,morphPoint);
        gRDMs = [{-RDMs(:,:,1)} {-RDMs(:,:,2:3)}]; % {-RDMs(:,:,3)}];
        [r2penalty tr2 shuffle_r2penalty shuffle_tr2 pop_r2penalty pop_tr2] = classFactorLoading( ...
            matchedSim.pearson,gRDMs,nchoosek(minSessions,2));
        across.class.d3.r2penalty = [across.class.d3.r2penalty; r2penalty];
        across.class.d3.tr2 = [across.class.d3.tr2; tr2];
        across.class.d3.pop.r2penalty = [across.class.d3.pop.r2penalty; pop_r2penalty];
        across.class.d3.pop.tr2 = [across.class.d3.pop.tr2; pop_tr2];
        across.class.d3.shuffle_r2penalty = [across.class.d3.shuffle_r2penalty; shuffle_r2penalty];
        across.class.d3.shuffle_tr2 = [across.class.d3.shuffle_tr2; shuffle_tr2];
        plotgRDMs(gRDMs,[root '/RDMs_3D']);

        % Drift vs Context Analysis

        partSim = partitionRDMs(matchedSim.pearson,envs);
        across.partSim = cat(5,across.partSim,partSim);
        across.sim = cat(3,across.sim,matchedSim.pearson);
        across.cellAnimalID = [across.cellAnimalID; mi.*ones(size(matchedSim.pearson,3),1)];

% %         % Context prediction stuff
        across.contextPrediction.rawPV.actual = cat(3,across.contextPrediction.rawPV.actual, ...
            decodeContext(matchedSim.pv,envs,6,6));
        across.contextPrediction.zeroedPV.actual = cat(3,across.contextPrediction.zeroedPV.actual, ...
            decodeContext(matchedSim.zeroedpv,envs,6,6));
    end
    
    mouseNames = [];
    for p = paths'
        i = slind(p{1});
        mouseNames = [mouseNames; {p{1}(i(end)+1:end-4)}];
    end
    
    root = ['Plots/BatchedAnalyses/Summary']; 
    set(0,'DefaultAxesColorOrder',plasma(7))';
    
    
    %%%%%%%%%%%%%%%% Context vs Drift %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    conXdrift(across.partSim,across.maps,across.envs, ...
        [root '/ContextVsDrift'],across.mfr);
    
    %%%%%%%%%%%%%%%%%%%%% GLM Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     for k = 1:3
%         colorMap = circshift([linspace(0,1,256)', zeros(256,2)],[0 k]);
%         colormap(colorMap);
%         colorbar;
%         saveFig(gcf,[root '/RDM_Colorbar_Examples_' num2str(k)],[{'tiff'} {'pdf'}]);
%     end
%     saveFig(gcf,[root '/RDM_Example_Inferno'],[{'tiff'} {'pdf'}]);

    tmp = across.class.d3.r2penalty; %bsxfun(@times,across.class.d2.tr2,across.class.d2.r2penalty);
    [maxTime maxTimeInds] = sort(tmp(:,1),'descend');
    maxTimeInds(isnan(maxTime)) = [];
    [maxContext maxContextInds] = sort(tmp(:,2),'descend');
    maxContextInds(isnan(maxContext)) = [];
    
    tmp = [];
    for mi = 1:5
        tmp = [tmp nansum(across.cellAnimalID(maxContextInds)==mi)];
    end
    fprintf(sprintf('\n\tn = %i',[tmp]));
    
    [across.class.d3.tr2(maxContextInds(1),:) across.class.d3.r2penalty(maxContextInds(1),:)]
    
%     plotStackMaps(across.maps(:,:,maxTimeInds(14),:));
    plotStackMaps(across.maps(:,:,maxTimeInds(7),:));
    saveFig(gcf,[root '/ExampleCells_Time'],[{'tiff'} {'pdf'}]);
    plotStackMaps(across.maps(:,:,maxContextInds(1),:));
    saveFig(gcf,[root '/ExampleCells_Context'],[{'tiff'} {'pdf'}]);
%     plotStackMaps(across.maps(:,:,maxShapeInds(7),:));
%     plotStackMaps(across.maps(:,:,maxShapeInds(10),:));
    
    pVal = plotGLM_2d((across.class.d2.r2penalty),...
        (across.class.d2.shuffle_r2penalty),[root '/GLM_Analysis/2D']);

    pVal = plotGLM_2d((across.class.d3.r2penalty),...
        (across.class.d3.shuffle_r2penalty),[root '/GLM_Analysis/3D'],[]);

%     plotGLM_3d(across,[root '/GLM_Analysis/3D']);
    %%%%%%%%%%%%%%%%%%%%% Context Prediction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     across = load('MapSVM');
%     plotDecoding(across.contextPrediction.justTracked.actual, ...
%         [root '/Decoding/Accuracy_TrackedSVM']);
    plotDecoding(across.contextPrediction.rawPV.actual, ...
        [root '/Decoding/Accuracy_RawPV']);
    plotDecoding(across.contextPrediction.zeroedPV.actual, ...
        [root '/Decoding/Accuracy_ZeroedPV']);
    
    %%%%%%%%%%%%%%%%%%%%% Place cell stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numPCSesh = nansum(across.isPlaceCell < 0.05,2);
    numSesh = nansum(~isnan(across.isPlaceCell),2);
    
    figure
    set(gcf,'position',[50 50 750 500])
    plot([0 33],[0 0],'linestyle','--','color','k')
    hold on
    plot([0 33],[1 1],'linestyle','--','color','k')
    nsims = 100;
    xl = repmat([1:32],[nsims 1]);
    tmp1 = bsxfun(@rdivide,cumsum(rand(nsims,32)<0.05,2),[1:1:32]);
    scatter(randn(length(xl(:)),1).*0.25+xl(:), ...
        randn(length(tmp1(:)),1).*0.025+tmp1(:),10,[0.5 0.5 0.5]);   
    scatter(randn(length(numSesh),1).*0.25+numSesh, ...
        randn(length(numSesh),1).*0.025+numPCSesh./numSesh,10,[0.2 0.2 0.9])    
    set(gca,'ylim',[-0.2 1.2],'xlim',[0 33],'xtick',[0:4:32])
    ylabel('Proportion sessions meeting place cell criteria')
    xlabel('Number of tracked sessions')
    saveFig(gcf,[root '/PlaceCellCriteria'],[{'tiff'} {'pdf'}])
    
    %%%%%%%%%%%%%%%%%%%%% Sampling Coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure
    set(gcf,'position',[50 50 550 200])
    h = plot(across.coverage','marker','o','markerfacecolor','auto');
    hold on
    set(gca,'xlim',[0 33],'xtick',[0:4:32],'ylim',[0 1])
    legend(h,mouseNames,'location','southeast','NumColumns',2)
    xlabel('Session Number')
    ylabel('Proportion sampled')
    saveFig(gcf,[root '/SamplingCoverage'],[{'tiff'} {'pdf'}])
    
    %%%%%%%%%%%%%%%%%%%%% MDS Quantification Plots %%%%%%%%%%%%%%%%%%%%%%
    
    figure
    set(gcf,'position',[50 50 150 200])
    h = plot(across.mds.driftAmount','marker','o','markerfacecolor','auto');
    hold on
    set(gca,'xlim',[0.5 5.5],'xtick',[1:5])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    saveFig(gcf,[root '/MDS_DriftAmount'],[{'tiff'} {'pdf'}])
    
    %%%%%%%%%%%%%%%%%%% Sequence Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure
    set(gcf,'position',[50 50 200 200])
    h = plot(across.seqAnalysis.transitionPoint'-3.5,'marker','o','markerfacecolor','auto');
    hold on
    set(gca,'xlim',[0.5 5.5],'ylim',[-2.5 2.5],'ytick',[-2.5:1:2.5],'xtick',[1:5])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    legend(h,mouseNames,'location','southwest')
    xlabel('Sequence')
    saveFig(gcf,[root '/SequenceAnalysis_TransitionPoint'],[{'tiff'} {'pdf'}])
    
    figure
    set(gcf,'position',[50 50 200 200])
    plot(across.seqAnalysis.maxDecorr','marker','o','markerfacecolor','auto')
    hold on
    set(gca,'xlim',[0.5 5.5],'ylim',[0 nanmax(get(gca,'ylim'))],'xtick',[1:5])
    xlabel('Sequence')
%     plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    saveFig(gcf,[root '/SequenceAnalysis_MaxDecorr'],[{'tiff'} {'pdf'}])
    
    [h pval ci tstat] = ttest(across.seqAnalysis.maxDecorr(:,1),across.seqAnalysis.maxDecorr(:,5))
end














































