function continuousMapAnalysis(paths,doPlot)
    
    clc
    close all
    drawnow
    
    doMDS = true;
    
    pause_thresh = 2;
    cell_thresh_for_burst = 0.1;
    doPreplay = false;
    doSubBins = false;
    subBinSize = 150;
    doVectorSim = false;
    doAnatomy = false;
    pvWindow = 60; % Frames for chunking vec for doVectorSim
    doCOM = false;
    if nargin < 2
        doPlot = false;
    end

    warning off all
    
    piece = [];
    spiece = [];
    amid = []; 
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
        blah = find(ismember(piece{end},'_'),2,'last');
        amid = [amid; {piece{end}(find(ismember(piece{end},'/'),1,'last')+1:blah(1)-1)}];
    end
    uamid = unique(amid);
    upiece = unique(piece);
    
    sigParams = [];
    envSize = [17 17];
    includeLags = inf;
    pThresh = 0.025;
    envLabel = [{'sq1'} {'sq2'} {'sq3'} ...
        {'g3'} {'g2'} {'g1'}];
    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        
        if ~isfield(s,'alignment')
            fprintf(['\n\t\t***** No full sequence continuous alignment for these sessions *****\n'])
            continue
        end
        
        doAl = help_getAlignmentID(s.alignment,length(isM),paths(isM));
        if isnan(doAl)
            fprintf(['\n\t\t***** No full sequence continuous alignment for these sessions *****\n'])
            continue
        end
        alignMap = s.alignment(doAl).alignmentMap;
        am = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        simVecs = repmat({[]},[1 length(sessions)]);
        aos = [];
        envs = [];
        tic
        
        slashInds = find(ismember(paths{1},'/'));
        root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
        slashInds = find(ismember(upiece{mi},'/'));
        root = [root '/' upiece{mi}(slashInds(end)+1:end) '/ContinuousSequenceAnalyses'];
        
        sequenceNum = str2num(upiece{mi}(find(ismember(upiece{mi},'_'),1,'last')+1:end));
        blah = find(ismember(upiece{mi},'_'),2,'last');
        mouseID = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:blah(1)-1);
        mouseID = find(ismember(uamid,mouseID));
        fprintf(['\t\tPreloading Data... '])
        for si = 1:length(sessions)
            if doAnatomy
                s = load(sessions{si},'processed','exclude','calcium');
            else
                s = load(sessions{si},'processed','exclude');
            end
            slashInds = find(ismember(sessions{si},'/'));
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            amfr{si} = nanmean(gT,2);
            v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            [m os] = mkTraceMaps(s.processed.p,gT,v>=pause_thresh,envSize);
            aos = cat(3,aos,os);
            am{si} = m;
            
            isPC{si} = s.processed.splithalf.wholemap_unmatched.p <= pThresh;
            if isfield(s.processed,'exclude')
                isPC{si} = isPC{si}(s.processed.exclude.SFPs,:);
            end
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
        end  
        durat = toc;
        fprintf([num2str(durat) ' s\n']);
        
        um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
        upc = nan([length(alignMap{1}(:,1)) length(sessions)]);
        for si = 1:length(sessions)
            isCPC = isPC{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
            tmp = false(length(alignMap{1}(:,si)),1);
            tmp(alignMap{1}(:,si)~=0) = isCPC;
            upc(:,si) = tmp; %%% store isPC

            um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));       
        end
        
        sim = repmat({[]},length(sessions));
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing pairwise map comparisons: '])
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
                fprintf(str);
                strLength = length(str);
                
                tmp1 = um(:,:,:,si);
                tmp2 = um(:,:,:,sj);
                
                isGood = upc(:,si)|upc(:,sj); %%% Include only place cells
                tmp1(:,:,~isGood) = nan;
                tmp2(:,:,~isGood) = nan;
                
                ex = permute(all(all(isnan(tmp1),1),2) | all(all(isnan(tmp2),1),2),[3 1 2]);
                tmp1(:,:,ex) = [];
                tmp2(:,:,ex) = [];
                
                %%% Choose best orientation
                
                vals = nan(1,4);
                for rot = 0:3
                    rtmp2 = imrotate(tmp2,rot.*90);
                    goodPixels = ~isnan(rtmp2(:,:,1))&~isnan(tmp1(:,:,1));
                    vals(rot+1) = corr(tmp1(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])), ...
                        rtmp2(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])));
                end
                [a bestRot] = nanmax(vals);
                
                tmp2 = imrotate(tmp2,(bestRot-1).*90);
                
                
                xc = pvxcorr3(tmp1,tmp2,[1 1 1],30);
                a = 2;
                b = 2;
                c = 2;
                ivals = xcorr3transform(tmp1,tmp2,[a-ceil(length(xc(:,1,1))./2) ...
                    b-ceil(length(xc(1,:,1))./2) c-ceil(length(xc(1,1,:))./2)]);
                sim{si,sj} = ivals;
            end
        end        
        
        help_showpairwise(sim);
        colormap jet;
        caxis([0 1]);
        saveFig(gcf,[root '/ConditionComparison_RDM'],[{'pdf'} {'tiff'}]);
        
%         cellfun(@length,sim)
%         cellfun(@nanmedian,sim)
        
        params = transitionPlot(sim,envs,doComps, ...
            [root '/ConditionComparison_Correlations_Cellwise']);
        
        if ~isempty(params)
            sigParams = [sigParams; params.sigmoidal_mad ...
                params.sigmoidal_intercept-3.5 params.sigmoidal_max_sep mouseID sequenceNum];
        end
        fprintf('\n')
        
        close all
    end
    
    figure
    set(gcf,'position',[50 50 600 300])
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,3);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    subplot(1,2,1)
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
%     plot(nanmean(toPlot)','marker','none','color','k','linewidth',2,'linestyle','-')
    set(gca,'xlim',[0 6])
%     h = mkGraph(toPlot);
    ylabel('Maximum decorrelation (Z(\Deltar))')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'ylim',[0 1.2])
    subplot(1,2,2)
    tmp = bsxfun(@minus,toPlot(:,2:end),toPlot(:,1));
    h = plot([zeros(length(toPlot(:,1)),1) tmp]', ...
        'marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
%     mkGraph(bsxfun(@minus,toPlot(:,2:end),toPlot(:,1)))
    set(gca,'xlim',[0 6])
%     h = mkGraph(toPlot);
    ylabel('Maximum decorrelation (Z(\Deltar))')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'ylim',[-0.4 0.4])
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    saveFig(gcf,[root(1:nanmin(find(ismember(root,'/'),2,'last'))-1) ...
        '/Maximum_Separation'],[{'pdf'} {'tiff'}]);
    
    
    
    figure
    set(gcf,'position',[50 50 300 300])
    hold on
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,2);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
    set(gca,'xlim',[0 6])
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    ylabel('median(|r_{sq1} - r_{g1}|)')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'ylim',[-2.5 2.5],'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('Transition Point')
    xlabel('Morph Sequence Number')
    saveFig(gcf,[root(1:nanmin(find(ismember(root,'/'),2,'last'))-1) ...
        '/SigFitIntercept'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 50 600 300])
    hold on
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,1);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    subplot(1,2,1)
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
    plot(nanmean(toPlot)','marker','none','color','k','linewidth',2,'linestyle','-')
    set(gca,'xlim',[0 6])
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    ylabel('median(|r_{sq1} - r_{g1}|)')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('mean(abs(sig-lin))')
    xlabel('Morph Sequence Number')
    subplot(1,2,2)
    h = plot([zeros(length(toPlot(:,1)),1) bsxfun(@minus,toPlot(:,2:end),toPlot(:,1))]', ...
        'marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    set(gca,'xlim',[0 6])
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    ylabel('median(|r_{sq1} - r_{g1}|)')
    xlabel('Sequence Number')
%     legend(h(:,1),uamid','location','eastoutside')
    set(gca,'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('mean(abs(sig-lin))')
    xlabel('Morph Sequence Number')
    
    saveFig(gcf,[root(1:nanmin(find(ismember(root,'/'),2,'last'))-1) ...
        '/MeanAbsDiffBtwn_SigLin'],[{'pdf'} {'tiff'}]);
end




















