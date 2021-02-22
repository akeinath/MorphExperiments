function pairwiseMapAnalysis(paths,doPlot)
    
    animalCondition = [{'CAMCA130'} {'Light'}; {'AKCA127'} {'Light'}; ...
        {'AKCA125'} {'Dark'}; {'AKCA170'} {'Dark'}];

    close all
    drawnow

    cell_thresh_for_burst = 0.1;
    pause_thresh = -2;
    
    
    doGeoAdjust = false;
    doPreplay = false;
    doVectorSim = false;
    doSubBins = false;
    doCOM = false;
    
    subBinSize = 150;
    pvWindow = 60;
    pThresh = 1; %0.025; % 0.025
    
    if nargin < 2
        doPlot = false;
    end

    warning off all
    
    clc
    fprintf('\nComputing pairwise map analysis ')
    if doPlot
        fprintf('and plotting maps')
    end
    fprintf('\n')
    %%% Reliability constaint
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    labels = [];
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
    
    includeLags = inf;
    envLabel = [{'sq1'} {'sq2'} {'sq3'} ...
        {'g3'} {'g2'} {'g1'}];
    doComps = [{'sq1'} {'sq1'}; {'sq1'} {'g1'}; {'sq1'} {'sq2'}; {'g1'} {'sq2'}; ...
        {'sq1'} {'sq3'}; {'g1'} {'sq3'}; {'sq1'} {'g3'}; {'g1'} {'g3'}; ...
        {'sq1'} {'g2'}; {'g1'} {'g2'}; {'sq1'} {'g1'}; {'g1'} {'g1'}]; %%% Complete curve
    
    aVecSim = [];
    geoSigParams = [];
    sigParams = [];
    environSize = [17 17];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)});
        doAl = help_getAlignmentID(s.alignment,2,paths(isM));
        if isnan(doAl)
            fprintf(['\n\t\t***** No pairwise alignment for these sessions *****\n'])
            continue
        end
        alignMap = s.alignment(doAl(1)).alignmentMap;
%         scoreMap = s.alignment(doAl(1)).scores;

        am = repmat({[]},[1 length(sessions)]);
        ap = repmat({[]},[1 length(sessions)]);
        abp = repmat({[]},[1 length(sessions)]);
        at = repmat({[]},[1 length(sessions)]);
        
        bam = repmat({[]},[20.*60.*30./(subBinSize.*30) length(sessions)]);
        amfr = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        aos = [];
        envs = [];
        preVecs = repmat({[]},[1 length(sessions)]);
        simVecs = repmat({[]},[1 length(sessions)]);
        COMs = repmat({[]},[1 length(sessions)]);
        tic
        fprintf(['\t\tGenerating maps... '])
        sequenceNum = str2num(upiece{mi}(find(ismember(upiece{mi},'_'),1,'last')+1:end));
        blah = find(ismember(upiece{mi},'_'),2,'last');
        mouseID = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:blah(1)-1);
        mouseID = find(ismember(uamid,mouseID));
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','exclude');
            slashInds = find(ismember(sessions{si},'/'));
%             fprintf(['\t\tPreloading for minimum sampling:  ' sessions{si}(slashInds(end)+1:end-4) '\n'])
            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            forV = imfilter(s.processed.p,fspecial('gauss',[1 150],30),'same','replicate');
            
%             v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            v = [0 sqrt(nansum(diff(forV,[],2).^2,1))].*30;
            at{si} = gT;
            ap{si} = s.processed.p;
            [m os] = mkTraceMaps(s.processed.p,gT,[v>=pause_thresh],environSize);
            am{si} = m;
            if doSubBins
                for bi = 1:20.*60./subBinSize
                    inBin = false(1,length(s.processed.p(1,:)));
                    inBin(nanmin((bi-1).*subBinSize.*30+1:(bi).*subBinSize.*30,length(gT(1,:)))) = true;
                    [tm] = mkTraceMaps(s.processed.p,gT,inBin&[v>=pause_thresh],environSize);
                    bam{bi,si} = tm;
                end
            end
            amfr{si} = nanmean(gT(:,v>=pause_thresh),2);
            aos = cat(3,aos,os);
            isPC{si} = s.processed.splithalf.wholemap_unmatched.p <= pThresh;
            if isfield(s.processed,'exclude')
                isPC{si} = isPC{si}(s.processed.exclude.SFPs,:);
            end
%             isPC{si} = true(size(isPC{si}));
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
            if doPreplay
                [preVecs{si} pvals] = getBurstVectors(s.processed.p,gT, ...
                    pause_thresh,cell_thresh_for_burst,isPC{si});
            end
            
            if doVectorSim
                simVecs{si} = chunkVec(gT,pvWindow);
                abp{si} = chunkVec(ap{si},pvWindow);
            end
            
            if doCOM
                COMs{si} = getCOMs(m);
            end
        end  
        
        durat = toc;
        fprintf([num2str(durat) ' s\n']);

        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                if si == sj
                    alignMap{si,sj} = [[1:length(am{si}(1,1,:))]' [1:length(am{si}(1,1,:))]'];
                end                
                if isempty(alignMap{si,sj})
                    alignMap{si,sj} = alignMap{sj,si}(:,[2 1]);
%                     scoreMap{si,sj} = scoreMap{sj,si};
                end
                if isempty(alignMap{si,sj})
                    continue
                end
                alignMap{si,sj} = alignMap{si,sj}(all(alignMap{si,sj}~=0,2),:); % scoreMap{si,sj}>=-0.95
            end
        end

        
        tic
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing map comparisons...\n\t\t\t'])
        crossSim = nan(length(sessions));
        kCrossCOM = repmat({[]},length(sessions));
        kCrossSim = repmat({[]},length(sessions));
        kDeformSim = repmat({[]},length(sessions));
        kCrossMFR = repmat({[]},length(sessions));
        binSim = nan([length(sessions).*(20.*60./subBinSize) length(sessions).*(20.*60./subBinSize)]);
        for si = 1:length(sessions)
            for sj = si:length(sessions)
                
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2)+length(sessions))]);
                
                fprintf(str);
                strLength = length(str);
                
                tmp1 = am{si};
                tmp2 = am{sj};
                
                if isempty(alignMap{si,sj})
                    continue
                end
                
                isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                
                if doCOM
                    kCrossCOM{si,sj} = [COMs{si}(alignMap{si,sj}(isGood,1),:) ...
                        COMs{sj}(alignMap{si,sj}(isGood,2),:)];
                end
                
                a = amfr{si}(alignMap{si,sj}(isGood,1));
                b = amfr{sj}(alignMap{si,sj}(isGood,2));
                kCrossMFR{si,sj} = abs(a-b)./nanmax(a,b);
                
%                 if nansum(isGood) < 1
%                     
%                 end
                
                tmp1 = tmp1(:,:,alignMap{si,sj}(isGood,1));
                tmp2 = tmp2(:,:,alignMap{si,sj}(isGood,2));
                otmp2 = tmp2;
                
                vals = nan(1,4);
                for rot = 0:3
                    rtmp2 = imrotate(tmp2,rot.*90);
                    goodPixels = ~isnan(rtmp2(:,:,1))&~isnan(tmp1(:,:,1));
                    vals(rot+1) = corr(tmp1(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])), ...
                        rtmp2(repmat(goodPixels,[1 1 length(tmp1(1,1,:))])));
                end
                [a bestRot] = nanmax(vals);
                
                tmp2 = imrotate(tmp2,(bestRot-1).*90);
                
% % %                 if length(alignMap{si,sj}(:,1))./nanmax(length(am{si}(1,1,:)), ...
% % %                     length(am{sj}(1,1,:))) < 0.25
% % %                     continue
% % %                 end
                
% % %                 xc = pvxcorr3(tmp1,tmp2,[3 3 5],20);
% % % %                 xc = pvxcorr3(tmp1,tmp2,[1 1 1],20);
% % %                 [a b c] = ind2sub(size(xc),nanmedian(find(xc==nanmax(xc(:)))));
% % %                 ivals = xcorr3transform(tmp1,tmp2,[a-ceil(length(xc(:,1,1))./2) ...
% % %                     b-ceil(length(xc(1,:,1))./2) c-ceil(length(xc(1,1,:))./2)]);
% % %                 
% % %                 crossSim(si,sj) = nanmax(xc(:));


                %%% Don't do pvxcorr correction
                xc = pvxcorr3(tmp1,tmp2,[1 1 1],30);
                a = 2;
                b = 2;
                c = 2;
                ivals = xcorr3transform(tmp1,tmp2,[a-ceil(length(xc(:,1,1))./2) ...
                    b-ceil(length(xc(1,:,1))./2) c-ceil(length(xc(1,1,:))./2)]);
                crossSim(si,sj) = xc(2,2,2);
                kCrossSim{si,sj} = ivals;
                
                if doGeoAdjust
                    
%                     kDeformSim{si,sj} = geoPathFits(ap{si},ap{sj}, ...
%                         at{si}(alignMap{si,sj}(isGood,1),:), ...
%                         at{sj}(alignMap{si,sj}(isGood,2),:));
                    
                    kDeformSim{si,sj} = fitDeformation(tmp1,tmp2); % OLD IMAGE-BASED RESCALING
                end                
                
                if doSubBins
                    for bi = 1:(20.*60./subBinSize)
                        for bj = 1:(20.*60./subBinSize)
                            tmp1 = bam{bi,si};
                            tmp2 = bam{bj,sj};
                            tmp1 = tmp1(:,:,alignMap{si,sj}(isGood,1));
                            tmp2 = tmp2(:,:,alignMap{si,sj}(isGood,2));
                            xc = pvxcorr3(tmp1,tmp2,[1 1 1],20);
                            binSim((si-1).*(20.*60./subBinSize)+bi, ...
                                (sj-1).*(20.*60./subBinSize)+bj) = xc(2,2,2);
                        end
                    end
                end
                
                if doPlot
                    toPlot = [tmp1 nan(length(tmp1(:,1,1)),4,nansum(isGood)) tmp2];

                    doK = [8 4];

                    for part = 0:floor(length(toPlot(1,1,:,1))/prod(doK))

                        figure(1)
                        set(gcf,'position',[50 50 900 1350])
                        for k = 1:prod(doK)
                            if part.*prod(doK)+k > length(toPlot(1,1,:))
                                break
                            end

                            subplot(doK(1),doK(2),k)
                            imagesc(toPlot(:,:,part.*prod(doK)+k))
                            colormap jet
    %                         caxis([0 nanmax(nanmax(toPlot(:,:,part.*prod(doK)+k)))])
                            alpha(double(~isnan(toPlot(:,:,part.*prod(doK)+k))))
                            axis equal
                            axis off    
                        end
                        slashInds1 = find(ismember(sessions{si},'/'));
                        slashInds2 = find(ismember(sessions{sj},'/'));
                        outP = ['Plots/PairwiseAlignedCellMaps_FromMapAnalysis/' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '/' ...
                            [sessions{si}(slashInds1(end)+1:end-4) '_vs_' sessions{sj}(slashInds2(end)+1:end-4)] ...
                            '_Part_' num2str(part)];
                        saveFig(gcf,outP,[{'tiff'} {'pdf'}])
                        close all
                        drawnow
                    end
                end
            end
        end
        durat = toc;
        fprintf(['\t' num2str(durat) ' s\n']);
        
        if doPreplay
            tic
            fprintf(['\t\tComparing burst vectors... '])
            preplaySim = repmat({[]},[length(sessions) length(sessions)]);
            for si = 1:length(sessions)
                vec = preVecs{si};
                for sj = 1:length(sessions)
                    if isempty(alignMap{si,sj}) 
                        continue
                    end
                    isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                        
                    tmp1 = vec(alignMap{si,sj}(isGood,1),:);
                    tmp2 = am{sj}(:,:,alignMap{si,sj}(isGood,2));
                    cv2 = permute(tmp2,[3 1 2]);
                    cv2 = reshape(cv2,[length(cv2(:,1,1)) prod(size(tmp2(:,:,1)))]);

                    rmat = corr(tmp1,cv2);
                    preplaySim{si,sj} = nanmax(rmat,[],2);
                end
            end
            durat = toc;
            fprintf([num2str(durat) ' s\n']);
        end
        
        if doVectorSim
            tic
            fprintf(['\t\tComparing vector similarity... '])
            vectorSim = repmat({[]},[length(sessions) length(sessions)]);
            
            isGoodCount = nan(length(sessions));
            for si = 1:length(sessions)
                for sj = si+1:length(sessions)
                    if isempty(alignMap{si,sj}) 
                        continue
                    end
                    isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                    isGood = true(length(alignMap{si,sj}(:,1)),1);    
                    isGoodCount(si,sj) = nansum(isGood);
                end
            end
            minIsGood = nanmin(isGoodCount(:));
            
            for si = 1:length(sessions)
%                 p1 = floor(abp{si}./2.5)+1;
                for sj = si+1:length(sessions)
                    if isempty(alignMap{si,sj}) 
                        continue
                    end
                    isGood = isPC{si}(alignMap{si,sj}(:,1)) | isPC{sj}(alignMap{si,sj}(:,2));
                    isGood = true(length(alignMap{si,sj}(:,1)),1);    
                    
                    newIsGood = false(size(isGood));
                    gInds = find(isGood);
                    gInds = gInds(randperm(length(gInds)));
                    newIsGood(gInds(1:minIsGood)) = true;
                    isGood = newIsGood;
                    
                    
                    tmp1 = simVecs{si}(alignMap{si,sj}(isGood,1),:);
                    tmp2 = simVecs{sj}(alignMap{si,sj}(isGood,2),:);
                    xc = corr(tmp1,tmp2);

                    
                    vectorSim{si,sj} = nanmax(xc,[],2);
                    vectorSim{sj,si} = nanmax(xc,[],1)';
                end
                vectorSim{si,si} = nan(nanmax(cellfun(@length,vectorSim(si,:))),1);
            end
            bmMat = nan(length(vectorSim));
            for i = 1:length(vectorSim)
                tmp = cat(2,vectorSim{i,:});
                [a b] = nanmax(tmp,[],2);
                for j = 1:length(vectorSim)
                    bmMat(j,i) = nanmean(b==j);
                end
            end
            [a condGroup] = ismember(envs,envLabel);
            tmp = repmat({[]},[nanmax(condGroup) nanmax(condGroup)]);
            for i = 1:length(sessions)
                for j = 1:nanmax(condGroup)
                    tmp{condGroup(i),j} = [tmp{condGroup(i),j}; ...
                        bmMat(condGroup==j,i)];
                end
            end
            vs = cellfun(@nanmean,tmp);
            close all
            figure
            imagesc(vs)
            colorbar
            drawnow
            durat = toc;
            fprintf([num2str(durat) ' s\n']);
        end
        
%         imagesc(crossSim)
%         colorbar

        slashInds = find(ismember(paths{1},'/'));
        root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
        slashInds = find(ismember(upiece{mi},'/'));
        root = [root '/' upiece{mi}(slashInds(end)+1:end)];

%         if length(unique(envs))<3
        help_showpairwise(kCrossSim);
        colormap jet;
        caxis([0 1]);
        cellfun(@nanmedian,kCrossSim)
%         cellfun(@nanmedian,kDeformSim)
        cellfun(@length,kCrossSim)
%             continue
%         end

%         help_showpairwise(kDeformSim);
%         colormap jet;
%         caxis([0 1]);


        if length(unique(envs))<3
            continue
        end
        
        params = transitionPlot(kCrossSim,envs,doComps, ...
            [root '/ConditionComparison_Correlations_Cellwise']);

        
        transitionPlot(crossSim,envs,doComps, ...
                [root '/ConditionComparison_PV']);
%             
%         transitionPlot(kCrossMFR,envs,doComps, ...
%                 [root '/ConditionComparison_MFR']);
            
%         mds2D(kCrossSim,envs,envLabel,[root '/Full_Cellwise_MDS_2D']);
            
        if length(unique(envs))<6
            continue
        end
            
        if ~isempty(params)
            sigParams = [sigParams; nanmean(params.sigmoidal_params(:,2)) ...
                params.poly_intercept-3.5 params.attractor_strength mouseID sequenceNum];
        end
        
        [a condGroup] = ismember(envs,envLabel);

%         help_showpairwise(groupMat(kCrossSim,condGroup));
%         saveFig(gcf,[root '/Sessionwise_Correlogram'],[{'pdf'} {'tiff'}]);
%         
%         help_showpairwise(groupMat(kCrossMFR,condGroup))
%         saveFig(gcf,[root '/Sessionwise_MFR_changes'],[{'pdf'} {'tiff'}]);

        if doSubBins
            mds2D_binned(binSim,envs,envLabel, ...
                (20.*60./subBinSize),[root '/Full_Cellwise_MDS_2D_TimeBinned']);
        end
                
        if length(kCrossSim(:,1,1))>8
%             mds3D(kCrossSim,envs,envLabel,[root '/Full_Cellwise_MDS_3D.gif']);
%             mds3D(kCrossMFR,envs,envLabel,[root '/Full_Cellwise_RateDiffs_MDS_3D.gif']);
        end
        
        if doPreplay
%             
%             for gi = 1:5
%                 goodLag = false(size(preplaySim)); 
%                 goodLag((gi-1).*6+1:(gi).*6,(gi-1).*6+1:(gi).*6) = true;
                doComps2 = [{'sq1'} {'sq2'} {'sq3'} {'g3'} {'g2'} {'g1'}];
                toPlot = repmat({[]},[2 length(doComps2)]);
                for i = 1:length(doComps2)
                    a = cat(1,preplaySim{logical(bsxfun(@times, ...
                        ismember(envs,doComps2(i)),ismember(envs,{'g1'})'))});
                    b = cat(1,preplaySim{logical(bsxfun(@times, ...
                        ismember(envs,doComps2(i)),ismember(envs,{'sq1'})'))});
                    toPlot{1,i} = a;
                    toPlot{2,i} = b;
                end
%                 cellfun(@nanmean,toPlot)
%             end
            
            
            figure
            set(gcf,'position',[50 50 250 350])
            mkGraph(toPlot)
            set(gca,'ylim',[-2 2])
            hold on
            plot(get(gca,'xlim'),[0 0],'linewidth',1,'linestyle','--','color','k')
            xlabel('Condition Comparison')
            
            saveFig(gcf,[root '/Burst_Similarity'],[{'pdf'} {'tiff'}]);
        end
        
        if doGeoAdjust
%             help_showpairwise(groupMat(kDeformSim,condGroup))
%             saveFig(gcf,[root '/GeoAdjusted/Sessionwise'],[{'pdf'} {'tiff'}]);
            
            params = transitionPlot(kDeformSim,envs,doComps, ...
                [root '/GeoAdjusted/ConditionComparison_Correlations_Cellwise']);
            
            if ~isempty(params)
                geoSigParams = [geoSigParams; nanmean(params.sigmoidal_params(:,2)) ...
                    params.poly_intercept-3.5 params.attractor_strength mouseID sequenceNum];
            end
            
            %%% Multidimensional 3D video
            if length(kDeformSim(:,1,1))>8
                mds2D(kDeformSim,envs,envLabel,[root '/GeoAdjusted/Full_Cellwise_MDS_2D']);
%                 mds3D(kDeformSim,envs,envLabel,[root '/GeoAdjusted/Full_Cellwise_MDS_3D.gif']);
            end
        end
        
        if doCOM
            comVecPlot(groupMat(kCrossCOM,condGroup,false),[root '/COM_Vectors']);
        end
        
        close all
        drawnow
    end
    
    figure
    set(gcf,'position',[50 50 500 500])
    toPlot = [];
    for i = 1:length(sigParams(:,end-1))
        toPlot{sigParams(i,end-1),sigParams(i,end)} = sigParams(i,3);
    end
    toPlot(cellfun(@isempty,toPlot)) = {nan};
    toPlot = cell2mat(toPlot);
    
    h = plot(toPlot','marker','o','markerfacecolor','auto','linewidth',2,'markersize',10);
    hold on
    plot(nanmean(toPlot)','marker','none','color','k','linewidth',2,'linestyle','-')
    set(gca,'xlim',[0 6])
%     h = mkGraph(toPlot);
    ylabel('Median decorrelation (\Deltar)')
    xlabel('Sequence Number')
    legend(h(:,1),uamid','location','northwest')
    set(gca,'ylim',[-0.05 0.7])
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--','linewidth',1)
    saveFig(gcf,[root(1:find(ismember(root,'/'),1,'last')-1) ...
        '/MedianMapDecorrelation'],[{'pdf'} {'tiff'}]);
    
    figure
    set(gcf,'position',[50 50 500 500])
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
    legend(h(:,1),uamid','location','northwest')
    set(gca,'ylim',[-2.5 2.5],'xtick',[1:1:nanmax(sigParams(:,5))])
    ylabel('Transition Point')
    xlabel('Morph Sequence Number')
    saveFig(gcf,[root(1:find(ismember(root,'/'),1,'last')-1) ...
        '/SigFitIntercept'],[{'pdf'} {'tiff'}]);

    if exist('geoSigParams')~=0 && ~isempty(geoSigParams) && length(geoSigParams(:,1)) > 1
        
        figure
        set(gcf,'position',[50 50 250 250])
        h(1) = plot(sigParams(:,end),sigParams(:,2),'color','k','linewidth',1.5);
        hold on
        if doGeoAdjust
            toPlot = [];
            for i = 1:length(geoSigParams(:,end-1))
                toPlot{geoSigParams(i,end-1),geoSigParams(i,end)} = geoSigParams(i,2);
            end
            toPlot = cellfun(@nanmean,toPlot);
            plot(repmat(1:length(toPlot(1,:)),[length(toPlot(:,1)) 1])', ...
                toPlot','color',[0.3 0.9 0.3],'linewidth',1.5,'marker','o',...
                'markerfacecolor',[1 1 1],'markersize',8);
            h(2) = plot(geoSigParams(:,end),geoSigParams(:,2),'color',[0.3 0.9 0.3],'linewidth',1.5);
            plot([0 6],[0 0],'color','k','linestyle','--')
        end
        set(gca,'ylim',[-2.5 2.5],'xtick',[1:1:nanmax(sigParams(:,5))])
        ylabel('Sigmoidal fit intercept')
        xlabel('Morph Sequence Number')
        if doGeoAdjust
            legend(h,[{'Raw'} {'Geometry-adjusted'}],'location','southwest')
        end
        saveFig(gcf,[root '/GeoAdjusted/SigFitIntercept'],[{'pdf'} {'tiff'}]);
        
        figure
        set(gcf,'position',[50 50 300 300])
        toPlot = [];
        for i = 1:length(geoSigParams(:,end-1))
            toPlot{geoSigParams(i,end-1),geoSigParams(i,end)} = geoSigParams(i,3);
        end
        h = mkGraph(toPlot);
        ylabel('median(|r_{sq1} - r_{g1}|)')
        xlabel('Sequence Number')
        legend(h(:,1),uamid','location','northwest')
        set(gca,'ylim',[-0.05 0.5])
        hold on
        plot([0 6],[0 0],'color','k','linestyle','--')
        saveFig(gcf,[root(1:find(ismember(root,'/'),1,'last')-1) ...
            '/GeoAdjusted/MedianMapDecorrelation'],[{'pdf'} {'tiff'}]);
% % % 
% % %         tmp = cellfun(@nanmean,toPlot);
% % %         tmp = diff(tmp,[],2);
% % %         [a b] = ismember(uamid,animalCondition(:,1)');
% % %         td = tmp(b(ismember(animalCondition(:,2),{'Dark'})),:);
% % %         tl = tmp(b(ismember(animalCondition(:,2),{'Light'})),:);
% % %         mkGraph([{td(:)} {tl(:)}])
        
%         figure
%         set(gcf,'position',[50 50 250 350])
%         h(1) = plot(sigParams(:,1),'color','k','linewidth',1.5);
%         hold on
%         if doGeoAdjust
%             h(2) =plot(geoSigParams(:,1),'color',[0.3 0.9 0.3],'linewidth',1.5);
%             plot([0 6],[0 0],'color','k','linestyle','--')
%         end
%         set(gca,'ylim',[0 10],'xtick',[1:1:length(sigParams(:,1))])
%         ylabel('Sigmoidal Steepness')
%         xlabel('Morph Sequence Number')
%         if doGeoAdjust
%             legend(h,[{'Raw'} {'Geometry-adjusted'}],'location','southoutside')
%         end
%         saveFig(gcf,[root '/SigFitSlope'],[{'pdf'} {'tiff'}]);
    end
end


















