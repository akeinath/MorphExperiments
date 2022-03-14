function batchAnalyses(paths)

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
        load(paths{mi},'um','umfr','upfr','ushc','envs','uIsPC','SFPs','aSamp');
        slashInds = find(ismember(paths{mi},'/'));
        root = ['Plots/BatchedAnalyses/' paths{mi}(slashInds(3)+1:end-4)];
        
        relativeCoverage = nan(1,length(aSamp(1,1,:)));
        for i = 1:32
            relativeCoverage(i) = nansum(nansum(aSamp(:,:,i)>=0.5,1),2) ./ ...
                nansum(nansum(any(aSamp(:,:,ismember(envs,envs(i)))>0,3),1),2);
        end
        across.coverage = [across.coverage; relativeCoverage];
        
% % % % % % % % % %         sim = getPairwiseMapSim(um,'pearson');
% % % % % % % % % %         sim.pfr = normCrop(upfr);
% % % % % % % % % %         sim.mfr = normCrop(umfr);
% % % % % % % % % %         sim.shc = normCrop(ushc);
% % % % % % % % % %         across.mfr = [across.mfr; umfr];
% % % % % % % % % %         
% % % % % % % % % % %         % SFPs similarity versus representational similarity
% % % % % % % % % % %         [sim.SFPs vssim counts] = sfpSimVS(SFPs,sim.pearson,minimum_pairwise_comparisons,[]);
% % % % % % % % % % %         [sim.SFPs vssim counts] = sfpSimVS(SFPs,sim.pearson,minimum_pairwise_comparisons,[root '/SFPsim_vs_PearsonRDM']);
% % % % % % % % % %         
% % % % % % % % % % %         [a b c] = mds2D(sim.pearson,envs,envLabel,[root '/MDS_2D']);
% % % % % % % % % %         [a b c] = mds2D(sim.pearson,envs,envLabel,[root '/MDS_2D']);
% % % % % % % % % %         across.mdsAngles = [across.mdsAngles; a];
% % % % % % % % % %         across.mdsDriftVariability = [across.mdsDriftVariability; b];
% % % % % % % % % %         across.mdsDriftAmount = cat(3,across.mdsDriftAmount,c); 
% % % % % % % % % %         [transitionPoint maxDecorr] = sim2seq(sim.pearson,uIsPC<0.05,envs,doComps,root);
% % % % % % % % % %         across.seqAnalysis.maxDecorr = [across.seqAnalysis.maxDecorr; maxDecorr];
% % % % % % % % % %         across.seqAnalysis.transitionPoint = [across.seqAnalysis.transitionPoint; transitionPoint];
% % % % % % % % % %         
% % % % % % % % % % %         plotSeqMaps(um,envs,envLabel,SFPs,root);
% % % % % % % % % %         
% % % % % % % % % %         % mean similarity across same context for pairwise seq comparison
% % % % % % % % % %         [seqSim simXlag] = conSimInSeq(sim,envs,envLabel);
% % % % % % % % % %         across.simXseq = cat(3,across.simXseq,seqSim);
% % % % % % % % % %         
% % % % % % % % % %         [conSim conXlag] = conComparisonInSeq(sim,envs,envLabel);
% % % % % % % % % %         across.simXcon = cat(3,across.simXcon,conSim);
% % % % % % % % % %         
% % % % % % % % % %         across.RDMs.partitioned = cat(5,across.RDMs.partitioned,partitionRDMs(sim,envs,envLabel));
% % % % % % % % % %         across.RDMs.whole = cat(3,across.RDMs.whole,sim.pearson);
% % % % % % % % % %         across.RDMs.envs = cat(2,across.RDMs.envs,repmat(envs,[1 length(sim.pearson(1,1,:))]));
% % % % % % % % % %         
%         fullConCells = find(nansum(nansum(~isnan(conSim),1),2)==15);
%         doV = conSim(:,:,fullConCells);
%         [contextParams contextError] = fitContextPattern(doV,1);


% % % % % %         doMFR = nanmedian(umfr,2);
% % % % % % 
% % % % % %         tmp = partitionRDMs(sim,envs,envLabel);
% % % % % %         conStack = [];
% % % % % %         for i = 1:5
% % % % % %             conStack = cat(3,conStack,permute(tmp(:,:,i,i,:),[1 2 3 5 4]));
% % % % % %         end
% % % % % %         nConComps = permute(nansum(nansum(~isnan(conStack),1),2),[4 3 1 2]);
% % % % % %         isFullCon = nConComps==nanmax(nConComps(:));
% % % % % %         [x y] = ind2sub(size(nConComps),find(isFullCon));
% % % % % % 
% % % % % %         doV = nan(6,6,length(x));
% % % % % %         for k = 1:length(x)
% % % % % %             doV(:,:,k) = tmp(:,:,y(k),y(k),x(k));
% % % % % %         end
% % % % % %         
% % % % % %         [contextParams contextError] = fitContextPattern(doV,1);
% % % % % %         
% % % % % %         conMAE = nan(size(nConComps));
% % % % % %         for k = 1:length(x)
% % % % % %             conMAE(x(k),y(k)) = contextError(k);
% % % % % %         end
% % %         
% % %         fullConCells = find(~isnan(nanmean(conMAE,2)));
% % %         error_rank = normRank(nanmean(conMAE(fullConCells,:),2),'ascend')./length(fullConCells);
% % %         
% % % %         error_rank = normRank(doMFR(~isnan(nanmean(conMAE,2)),:));
% % %         
% % %         vals = [];
% % %         valSim = [];
% % %         inc = (1./3);
% % %         for thresh = inc:inc:1
% % %             [a b] = contextSVM(um(:,:,fullConCells(error_rank>thresh-inc & error_rank<=thresh),:),envs,6);
% % %             vals = [vals; a];
% % %             valSim = [valSim; b];
% % %         end
% % %         figure
% % %         set(gcf,'position',[50 50 400 300])
% % %         doC = v2rgb(1:length(vals(:,1))).*0.5;
% % %         for i = 1:length(vals(:,1))
% % %             plot(vals(i,:)','color',doC(i,:))
% % %             hold on
% % %             set(gca,'ylim',[0 1])
% % %         end
% % %         plot(get(gca,'xlim'),[0.5 0.5],'linestyle','--','color',[0.5 0.5 0.5])
% % %         drawnow
% % %         
% % %         across.contextPrediction.contextErrorDivisions  = ...
% % %             cat(3,across.contextPrediction.contextErrorDivisions,vals);
% % %         
% % %         across.contextPrediction.contextErrorDivisions_similarity = ...
% % %             cat(3,across.contextPrediction.contextErrorDivisions_similarity,valSim);
        
%         clusteredContext = clustParams(contextParams,[],'mahalanobis');
%         plotClustProps(clusteredContext,doV,[root '/Clustering_Stability/']);
%         save('Clustered_Context','clusteredContext');
%         getEmbeddingExamples(clusteredContext,doV,[root '/Clustering_Stability/']);        
        
        
        
% % % % % %         [numSDrift mfrRankDrift conRankDrift] = quantPopDrift(um,umfr,conMAE);
% % % % % %         
% % % % % %         across.drift.numSPart = cat(3,across.drift.numSPart,numSDrift);
% % % % % %         across.drift.mfrRankPart = cat(3,across.drift.mfrRankPart,mfrRankDrift);
% % % % % %         across.drift.conRankPart = cat(3,across.drift.conRankPart,conRankDrift);
% % % % % % 
% % % % % %         [a b] = contextSVM(um,envs,6);
% % % % % %         across.contextPrediction.actual = cat(1,across.contextPrediction.actual,a);
% % % % % %         across.contextPrediction.actual_similarity = ...
% % % % % %             cat(1,across.contextPrediction.actual_similarity,b);
% % % % % %         
% % % % % %         [a b] = environmentSVM(um,envs,6);
% % % % % %         across.environmentPrediction.actual = cat(1,across.environmentPrediction.actual,a);
% % % % % %         across.environmentPrediction.actual_similarity = ...
% % % % % %             cat(1,across.environmentPrediction.actual_similarity,b);
        
        
        
        
%         [excludedCellSim cell_ids embedding] = c2cSim(sim.pearson,minimum_pairwise_comparisons);
        
        
%         [morphPoint] = sim2seq(sim.pearson(:,:,context_clustering.inds{5}),[],envs,doComps,root);
        
        % Time vs Context coding factor analysis
% % % % % %         morphPoint = ceil(transitionPoint);
% % % % % %         morphPoint(end+1) = morphPoint(end);
% % % % % %         RDMs = getRDMs(envs,envLabel,morphPoint);
% % % % % %         
% % % % % %         useNSessions = 12;
% % % % % %         
% % % % % %         for i = 1:32-useNSessions
% % % % % %             useSessions = i:i+useNSessions;
% % % % % %             minimum_pairwise_comparisons = 50;
% % % % % %             [r2penalty tr2 shuffle_r2penalty shuffle_tr2] = classFactorLoading( ...
% % % % % %                 sim.pearson(useSessions,useSessions,:),[{-RDMs(useSessions,useSessions,1)} ...
% % % % % %                 {-RDMs(useSessions,useSessions,2:3)}], ...
% % % % % %                 minimum_pairwise_comparisons);
% % % % % % 
% % % % % %             tr2penalty = bsxfun(@times,r2penalty,tr2);
% % % % % %             str2penalty = bsxfun(@times,shuffle_r2penalty,shuffle_tr2);
% % % % % %         end
% % % % % %         
% % % % % %         across.class.d2.r2penalty = [across.class.d2.r2penalty; r2penalty];
% % % % % %         across.class.d2.tr2 = [across.class.d2.tr2; tr2];
% % % % % %         across.class.d2.shuffle_r2penalty = [across.class.d2.shuffle_r2penalty; shuffle_r2penalty];
% % % % % %         across.class.d2.shuffle_tr2 = [across.class.d2.shuffle_tr2; shuffle_tr2];
% % % % % %         
% % % % % %         [r2penalty tr2 shuffle_r2penalty shuffle_tr2] = classFactorLoading( ...
% % % % % %             sim.pearson(useSessions,useSessions,:),-RDMs,minimum_pairwise_comparisons);
% % % % % %         across.class.d3.r2penalty = [across.class.d3.r2penalty; r2penalty];
% % % % % %         across.class.d3.tr2 = [across.class.d3.tr2; tr2];
% % % % % %         across.class.d3.shuffle_r2penalty = [across.class.d3.shuffle_r2penalty; shuffle_r2penalty];
% % % % % %         across.class.d3.shuffle_tr2 = [across.class.d3.shuffle_tr2; shuffle_tr2];

    end
        
    mouseNames = [];
    for p = paths'
        i = slind(p{1});
        mouseNames = [mouseNames; {p{1}(i(end)+1:end-4)}];
    end
    
    root = ['Plots/BatchedAnalyses/Summary']; 
    
    figure(1)
    set(gcf,'position',[50 50 350 250])
    h = plot(across.coverage');
    hold on
    set(gca,'xlim',[0 33],'ylim',[0 1])
    legend(h,mouseNames,'location','southwest')
    xlabel('Session (day)')
    ylabel('Proportion of environment visited')
    saveFig(gcf,[root '/EnvironmentCoverage'],[{'tiff'} {'pdf'}])
    
    figure
    set(gcf,'position',[50 50 300 200])
    h = plot(across.seqAnalysis.transitionPoint'-3.5,'marker','o');
    hold on
    set(gca,'xlim',[0.5 5.5],'ylim',[-2.5 2.5],'ytick',[-2.5:1:2.5])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    legend(h,mouseNames,'location','southwest')
    saveFig(gcf,[root '/SequenceAnalysis_TransitionPoint'],[{'tiff'} {'pdf'}])
    
    figure
    set(gcf,'position',[50 50 300 200])
    plot(across.seqAnalysis.maxDecorr','marker','o')
    hold on
    set(gca,'xlim',[0.5 5.5],'ylim',[0 nanmax(get(gca,'ylim'))])
%     plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    saveFig(gcf,[root '/SequenceAnalysis_MaxDecorr'],[{'tiff'} {'pdf'}])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxDays = 26;
    binned_vals = plotContextPrediction(across.contextPrediction.actual, ...
        across.contextPrediction.actual_similarity, ...
        [root '/Prediction_Analyses/'],maxDays);
    
    
    [a b c d] = ttest(cat(2,binned_vals{:}),1./2);
    fprintf('\n t(%i) = %0.2f, p = %.2e',[d.df; d.tstat; b./2])
    
    [binned_vals binned_sim_vals] = plotEnvironmentPrediction(across.environmentPrediction.actual, ...
        across.environmentPrediction.actual_similarity, ...
        [root '/Prediction_Analyses/'],maxDays);
    
    [a b c d] = ttest(cat(2,binned_vals{:}),1./6);
    fprintf('\n t(%i) = %0.2f, p = %.2e',[d.df; d.tstat; b./2])
    
    plotContextPrediction(permute(across.contextPrediction.contextErrorDivisions,[3 2 1]), ...
        permute(across.contextPrediction.contextErrorDivisions_similarity,[3 2 1]), ...
        [root '/Prediction_Analyses/'],maxDays)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plotStackedLine(across.drift.numSPart)
   
    figure(1)
    set(gcf,'position',[50 50 300 225])
    h = plotStackedLine(across.drift.mfrRankPart);
    set(gca,'ylim',[-0.1 0.7])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    legend(h,[{'1st'} {'2nd'} {'3rd'} {'4th'}])
    
    
    figure(1)
    set(gcf,'position',[50 50 300 225])
    h = plotStackedLine(across.drift.conRankPart);
    set(gca,'ylim',[-0.1 0.7])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    legend(h,[{'1st'} {'2nd'} {'3rd'} {'4th'}])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp = permute(across.contextPrediction.contextErrorDivisions,[3 2 1]);
    toPlot = [];
    for i = 1:length(tmp(1,1,:))
        toPlot = [toPlot {tmp(:,1:maxDays,i)}];
    end
    
    figure
    set(gcf,'position',[50 50 400 300])
    doC = (1-v2rgb(1:length(toPlot)+1)).*0.5 + 0.25;
    doC(1,:) = [];
    h = mkLine(toPlot,[],flipud(doC));
    set(gca,'ylim',[0 1])
    hold on
    plot(get(gca,'xlim'),[0.5 0.5],'linestyle','--','color',[0.5 0.5 0.5])
    legend(h(1,:),[{'1st'} {'2nd'} {'3rd'} {'4th'}],'location','southeast')
    saveFig(gcf,[root '/Prediction_Analyses/ContextGroupPrediction_XContextualFitError'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    
    quantRDM(across,[root '/RDM_Analyses'])
    
    seqSimWMatching(across,root);
    rdmClusterExplosion(across,[root,'/Clustering/']);
    
%     figure()
%     set(gcf,'position',[50 50 200 200])
%     h = cumHist([{across.class.d3.tr2} {across.class.d3.shuffle_tr2}],[0:0.01:1]);
%     [sstr2] = sort(across.class.d3.shuffle_tr2(~isnan(across.class.d3.shuffle_tr2)));
%     cutoff = sstr2(round(length(sstr2).*0.95));
%     
%     ar2 = bsxfun(@times,across.class.d2.r2penalty,across.class.d2.tr2);
%     sar2 = bsxfun(@times,across.class.d2.shuffle_r2penalty,across.class.d2.shuffle_tr2);
%     plot2DFactorLoadings(ar2,sar2);
%     
%     ar3 = bsxfun(@times,across.class.d3.r2penalty,across.class.d3.tr2);
%     sar3 = bsxfun(@times,across.class.d3.shuffle_r2penalty,across.class.d3.shuffle_tr2);
%     
%     doInclude = across.class.d3.tr2>cutoff;
%     figure
%     set(gcf,'position',[50 50 300 300])
%     plotTernary(across.class.d3.r2penalty(doInclude,:),[], ...
%         factor3dcm(across.class.d3.r2penalty(doInclude,:)), ...
%         [{'Drift'} {'Context Group'} {'Shape'}]);
end
















