function lagVals = getConVsStab_Matched(conMAE,ss,doMatch,root)

    conVsMFR = [];
    for i = 1:5
        tmp = nanmean(doMatch(:,(i-1).*6+1:(i).*6),2);
        conVsMFR = [conVsMFR; conMAE(~isnan(conMAE(:,i)),i) tmp(~isnan(conMAE(:,i)))];
    end

    
    figure
    set(gcf,'position',[50 50 250.*5 250])
    subplot(1,5,1)
    scatter(log10(conVsMFR(:,1)),log10(conVsMFR(:,2)),5)
    lsline
    axis equal
    axis square
    tmp = conVsMFR;
    [rvals pvals] = corr(tmp(~any(isnan(tmp),2),1),tmp(~any(isnan(tmp),2),2),'type','spearman');
    textX = get(gca,'xlim');
    textX = textX(2) - [textX(2)-textX(1)].*0.95;
    textY = get(gca,'ylim');
    textY = textY(2) - [textY(2)-textY(1)].*0.85;
    text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals pvals nansum(~any(isnan(tmp),2))]), ...
        'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','left');
    xlabel('Context fit error log10(MSE)')
    ylabel('log10(MFR)')
    set(gca,'ylim',[-4 -1]);
    tss = mat2lag(ss);
        
    for doLag = 1:4
        subplot(1,5,doLag+1)
        scatter(atanh(tss(:,doLag)),log10(nanmean(doMatch,2)),5)
        set(gca,'xlim',[-1 4])
        lsline
        axis equal
        axis square
        tmp = [tss(:,doLag) nanmean(doMatch,2)];
        [rvals pvals] = corr(tmp(~any(isnan(tmp),2),1),tmp(~any(isnan(tmp),2),2),'type','spearman');
        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.05;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.85;
        text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals pvals nansum(~any(isnan(tmp),2))]), ...
            'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','right');
        xlabel('Stability atanh(r)')
        ylabel('log10(MFR)')
        set(gca,'ylim',[-4 -1]);
    end
    saveFig(gcf,[root '/Context_Drift_MFR_' num2str(doLag)],[{'tiff'} {'pdf'} {'jpeg'}])
    

    tmp = [tss(:,doLag) nanmean(doMatch,2)];
    [rvals pvals] = corr(tmp(~any(isnan(tmp),2),1),tmp(~any(isnan(tmp),2),2),'type','spearman');
    sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals pvals nansum(~any(isnan(tmp),2))])

    figure
    set(gcf,'position',[50 50 250.*5 250])
    subplot(1,5,1)
    scatter(log10(conVsMFR(:,1)),log10(conVsMFR(:,2)),5)
    lsline
    axis equal
    axis square
    tmp = conVsMFR;
    [rvals pvals] = corr(tmp(~any(isnan(tmp),2),1),tmp(~any(isnan(tmp),2),2),'type','spearman');
    textX = get(gca,'xlim');
    textX = textX(2) - [textX(2)-textX(1)].*0.95;
    textY = get(gca,'ylim');
    textY = textY(2) - [textY(2)-textY(1)].*0.85;
    text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals pvals nansum(~any(isnan(tmp),2))]), ...
        'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','left');
    xlabel('Context fit error log10(MSE)')
    ylabel('log10(MFR)')
    set(gca,'ylim',[-4 -1]);
    tss = mat2lag(ss);

    for doLag = 1:4
        subplot(1,5,doLag+1)
        scatter(atanh(tss(:,doLag)),log10(nanmean(doMatch,2)),5)
        set(gca,'xlim',[-1 4])
        lsline
        axis equal
        axis square
        tmp = [tss(:,doLag) nanmean(doMatch,2)];
        [rvals pvals] = corr(tmp(~any(isnan(tmp),2),1),tmp(~any(isnan(tmp),2),2),'type','spearman');
        textX = get(gca,'xlim');
        textX = textX(2) - [textX(2)-textX(1)].*0.05;
        textY = get(gca,'ylim');
        textY = textY(2) - [textY(2)-textY(1)].*0.85;
        text(textX,textY,sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rvals pvals nansum(~any(isnan(tmp),2))]), ...
            'fontweight','normal','fontname','arial','fontsize',10,'horizontalalignment','right');
        xlabel('Stability atanh(r)')
        ylabel('log10(MFR)')
        set(gca,'ylim',[-4 -1]);
    end
    saveFig(gcf,[root '/Context_Drift_MFR_' num2str(doLag)],[{'tiff'} {'pdf'} {'jpeg'}])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numGroups = 4;
    lagVals = repmat({[]},[numGroups 4]);
    alv = repmat({[]},[1 4]);
    for i = 1:length(conMAE(1,:))
        doSim = [1:i-1 i+1:length(conMAE(1,:))];
        ts = ss(doSim,doSim,:);
        lag = abs(bsxfun(@minus,doSim,doSim'));
        
        tcm = repmat(permute(conMAE(:,i),[3 2 1]),[length(doSim) length(doSim)]);
        tlm = repmat(lag(:),[1 1 length(tcm(1,1,:))]);
        
        isGoodCon = find(~isnan(conMAE(:,i)));
        [a b] = sort(conMAE(isGoodCon,i));
        sortedCon = isGoodCon(b);
        
        
        
        groupInds = repmat({[]},[1 numGroups]);
        for gi = 1:numGroups
            groupInds{gi} = [sortedCon(floor((gi-1).*(length(sortedCon)./numGroups))+1: ...
                floor((gi).*(length(sortedCon)./numGroups)))];
        end
    
        tmpMatch = nanmean(doMatch(:,(i-1).*6+1:(i).*6),2);
        gm = repmat({[]},[1 numGroups]);
        for gi = 1:numGroups
            gm{gi} = tmpMatch(groupInds{gi});
        end
        
        [blah matchedInds] = compHist(gm, nanmin(tmpMatch):range(tmpMatch)./20:nanmax(tmpMatch),false);

        matchedGroupInds = [];
        for gi = 1:numGroups
            matchedGroupInds = [matchedGroupInds groupInds{gi}(matchedInds(:,gi))];
        end
        
        figure
        set(gcf,'position',[50 50 300 400])
        subplot(2,1,1)
        compHist(gm, nanmin(tmpMatch):range(tmpMatch)./20:nanmax(tmpMatch),true);
        set(gca,'ylim',[0 30])
        subplot(2,1,2)
        compHist([{tmpMatch(matchedGroupInds(:,1),1)} {tmpMatch(matchedGroupInds(:,2),1)} ...
            {tmpMatch(matchedGroupInds(:,3),1)} {tmpMatch(matchedGroupInds(:,4),1)}], ...
            nanmin(tmpMatch):range(tmpMatch)./20:nanmax(tmpMatch));
        set(gca,'ylim',[0 30])
        saveFig(gcf,[root '/MatchingFor_MFR_Example_' num2str(i)],[{'tiff'} {'pdf'} {'jpeg'}])
        close(gcf);
        
        tcm = conMAE(isGoodCon,i);
        tts = reshape(ts(:,:,isGoodCon),[16 1 length(isGoodCon)]);
        ttlm = tlm(:,:,isGoodCon);
       
        tlv = nan(length(isGoodCon),4);
        for j = 1:4
            tmp = reshape(tts(ttlm==j),[nansum(ttlm(:,1)==j) length(isGoodCon)]);
            tlv(:,j) = nanmean(tmp);
        end
        
        for j = 1:4
            alv{j} = [alv{j}; tcm(~isnan(tlv(:,j))) tlv(~isnan(tlv(:,j)),j)];
        end
        
        for gi = 1:length(matchedGroupInds(1,:))
            tts = reshape(ts(:,:,matchedGroupInds(:,gi)),[16 1 length(matchedGroupInds(:,gi))]);
            ttlm = tlm(:,:,matchedGroupInds(:,gi));
            
            tlv = nan(length(matchedGroupInds(:,gi)),4);
            for j = 1:4
                tmp = reshape(tts(ttlm==j),[nansum(ttlm(:,1)==j) length(matchedGroupInds(:,gi))]);
                tlv(:,j) = nanmean(tmp);
            end
            
            for j = 1:4
                lagVals{gi,j} = [lagVals{gi,j}; tlv(~isnan(tlv(:,j)),j)];
            end
        end
    end
    
    figure
    set(gcf,'position',[50 50 250.*4 250])
    for j = 1:4
        [rval pval] = corr(log10(alv{j}(:,1)),atanh(alv{j}(:,2)),'type','spearman');
        sprintf(['rho = %.3f\np = %.2e\nn = %.0f'],[rval pval length(alv{j}(:,1))])
        subplot(1,4,j)
        scatter(log10(alv{j}(:,1)),atanh(alv{j}(:,2)),5)
        axis square
        set(gca,'ylim',[-1 3],'xlim',[-5 0])
        lsline
        hold on
        plot(get(gca,'xlim'),[0 0],'linestyle','--')
    end
    saveFig(gcf,[root '/Context_Drift_GeneralCorrelation'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure
    set(gcf,'position',[50 50 400 250])
    mkBar(lagVals);
    hold on
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
    set(gca,'ylim',[-1 1])
    xlabel('Lag (Sequence)')
    ylabel('Correlation (r)')
    saveFig(gcf,[root '/Context_Drift_MatchingFor_MFR'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    
end





























