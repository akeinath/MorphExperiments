function [numSDrift mfrRankDrift conRankDrift] = quantPopDrift(um,umfr,conMAE)

    isTracked = permute(~all(all(isnan(um),1),2),[3 4 1 2]);
    numS = nansum(isTracked,2);    
    
    inc = 8;
    numSDrift = nan(32./inc,31);
    for minS = 0:inc:32-inc
        pvc = getPVC(um(:,:,numS>minS & numS<=minS+inc,:));
        [a b] = mat2lag(pvc);
        numSDrift(round(minS./inc) + 1,:) = a;
    end
    
    
    figure
    set(gcf,'position',[50 50 300 250]);
    doC = (1-v2rgb(1:length(numSDrift(:,1)))).*0.75;
    for minS = 1:1:length(numSDrift(:,1))
        plot(1:31,numSDrift(minS,:),'color',doC(minS,:));
        hold on
    end
    set(gca,'ylim',[-0.2 0.7])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mfrRank = normRank(nanmean(umfr,2))./length(umfr(:,1));
    inc = 0.25;
    mfrRankDrift = nan(1./inc,31);
    for minS = 0:inc:1-inc
        pvc = getPVC(um(:,:,mfrRank>minS & mfrRank<=minS+inc,:));
        [a b] = mat2lag(pvc);
        mfrRankDrift(round(minS./inc) + 1,:) = a;
    end
   
    figure
    set(gcf,'position',[50 50 300 250]);
    doC = (1-v2rgb(1:length(mfrRankDrift(:,1)))).*0.75;
    for minS = 1:1:length(mfrRankDrift(:,1))
        plot(1:31,mfrRankDrift(minS,:),'color',doC(minS,:));
        hold on
    end
    set(gca,'ylim',[-0.2 0.7])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mfrRank = normRank(nanmean(conMAE,2))./length(umfr(:,1));
    inc = 0.25;
    conRankDrift = nan(1./inc,31);
    for minS = 0:inc:1-inc
        pvc = getPVC(um(:,:,mfrRank>minS & mfrRank<=minS+inc,:));
        [a b] = mat2lag(pvc);
        conRankDrift(round(minS./inc) + 1,:) = a;
    end
   
    figure
    set(gcf,'position',[50 50 300 250]);
    doC = (1-v2rgb(1:length(conRankDrift(:,1)))).*0.75;
    for minS = 1:1:length(conRankDrift(:,1))
        plot(1:31,conRankDrift(minS,:),'color',doC(minS,:));
        hold on
    end
    set(gca,'ylim',[-0.2 0.7])
    plot(get(gca,'xlim'),[0 0],'linestyle','--','color','k')
end