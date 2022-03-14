function lagVals = getConVsStab_Matched2(conMAE,ss,doMatch,root)

    conVsMFR = [];
    for i = 1:5
        tmp = nanmean(doMatch(:,(i-1).*6+1:(i).*6),2);
        conVsMFR = [conVsMFR; conMAE(~isnan(conMAE(:,i)),i) tmp(~isnan(conMAE(:,i)))];
    end

    figure
    set(gcf,'position',[50 50 500 250])
    subplot(1,2,1)
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
    subplot(1,2,2)
    doLag = 1;
    scatter(atanh(tss(:,doLag)),log10(nanmean(doMatch,2)),5)
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
    saveFig(gcf,[root '/Context_Drift_MFR'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    numGroups = 4;
    lagVals = repmat({[]},[numGroups 4]);
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
        
% % %         Check that firing rate distributions are matched
%         [blah matchedInds] = compHist([{tmpMatch(matchedGroupInds{1},j)} {tmpMatch(matchedGroupInds{2},j)}],[0:0.0025:0.08]);
        

        for gi = 1:length(matchedGroupInds(1,:))
            tts = ts(:,:,matchedGroupInds(:,gi));
            ttlm = tlm(:,:,matchedGroupInds(:,gi));
            vals = [tts(:) ttlm(:)];
            
            goodVals = vals(~any(isnan(vals),2),:);
            
            for j = 1:4
                lagVals{gi,j} = [lagVals{gi,j}; goodVals(j==round(goodVals(:,2)),1)];
            end
        end
    end
    
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