function pred2Figs()
    paths = getFilePaths('MatlabData/Model_Grouped_SVMs_6Days','.mat');

    allSel = [];
    allAcc = [];
    allCoef = [];
    for p = paths'
        
        root = ['Plots/Model/' p{1}(find(ismember(p{1},'/'),1,'last')+1:end-4)];
        
        load(p{1},'tr2','linAcc','r2penalty');
        [sr2 best] = sort(r2penalty,2,'descend');
        isAttractor = best(:,1)==2;
        selectivity = [sr2(:,1) - nanmean(sr2(:,2:3),2)].*tr2(:,1);
%         selectivity = [r2penalty(:,2)-r2penalty(:,1)];
        
%         close all
%         plotWeightClass(r2penalty,tr2, ...
%             [root '/Ternary.gif'],false);
%         hist(selectivity,[-1:0.05:1])
%         drawnow
        
        allSel = [allSel; nanmean(selectivity)];
        
        allAcc = [allAcc; linAcc];
        
        
        try
            x = 1:80;
            data = linAcc(x);
            tbl = table(x', data');
            modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1)) + b(3);  
            beta0 = [1, 1, 1]; 
            mdl = fitnlm(tbl, modelfun, beta0);
            coefficients = mdl.Coefficients{:, 'Estimate'};
            fitData = modelfun(coefficients,x');
            allCoef = [allCoef; coefficients'];  
        catch
            allCoef = [allCoef; nan(1,3)];
        end

%         figure(1)
%         set(gcf,'position',[50 50 250 250])
%         plot(x,data,'color','k');
%         hold on
%         plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--');
%         plot(x,fitData,'color','r');
%         set(gca,'ylim',[0 1])
    end
%     scatter(allSel(:,1),allAcc(:,2))
    
    close all
    
    x = 1:80;
    low = allAcc(1:80,:);
    high = allAcc(81:end,:);
    figure()
    set(gcf,'position',[50 50 400 250])
    h = mkLine([{low(:,x)} {high(:,x)}],x);
    set(gca,'ylim',[0 1],'xlim',[nanmin(x)-1 nanmax(x)+1],'xtick',x(5:6:end),'xticklabel',x(5:6:end));
    hold on
%     h2 = patch([train(1) train(end) train(end) train(1)],[0 0 1 1],[0.5 0.5 0.5],...
%         'facecolor',[0.5 0.5 0.5],'edgecolor','none');
    plot(get(gca,'xlim'),[0.5 0.5],'color','k','linestyle','--')
    legend(h,[{'Lower selectivity'} {'Higher selectivity'}],'location','southeast')    
    pvals = nan(1,length(x));
    for k = 1:length(x)
        [a b c] = ranksum(low(:,k),high(:,k));
        pvals(k) = a;
    end
    text(x(pvals<0.05),0.9.*ones(1,nansum(pvals<0.05)),'*','fontname','arial', ...
        'fontweight','normal','fontsize',9,'horizontalalignment','center')
    xlabel('Day')
    ylabel('Prediction accuracy (%)')
    set(gca,'ytick',[0:0.25:1],'yticklabels',[0:25:100])
    drawnow
    outP = ['Plots/Model/Summary/SlidingWindowDecoding_6Days'];
    saveFig(gcf,outP,[{'tiff'} {'pdf'}])
    
    figure()
    set(gcf,'position',[50 50 200 250])
    h = mkBow([{allSel(1:80)} {allSel(81:end)}],[{'Lower'} {'Higher'}]);
    xlabel('Model group')
    ylabel('Selecitivity index')
    outP = ['Plots/Model/Summary/Grouped_SelectivityDiff'];
    saveFig(gcf,outP,[{'tiff'} {'pdf'}])
    
    endColors = [1 0.6 0.2; 0.2 0.6 1];
    rSel = [allSel-nanmin(allSel)]./[nanmax(allSel)-nanmin(allSel)];
    rSel = ((rSel-0.5).*1.5)+0.5;
    rSel(rSel<0) = 0;
    rSel(rSel>1) = 1;
    doColor = bsxfun(@plus,endColors(1,:),bsxfun(@times,rSel,diff(endColors,[],1)));

    figure
    set(gcf,'position',[50 50 900 250]);
    doCoef = ~isnan(allCoef(:,1));
    for i = 1:3
        subplot(1,3,i)
        scatter(allSel(doCoef),allCoef(doCoef,i),15,doColor(doCoef,:), ...
            'filled');
        lsline
        addCorr(allSel(doCoef),allCoef(doCoef,i))
        xlabel('Selecitivity index')
        ylabel(['Model coef ' num2str(i)])
    end
    outP = ['Plots/Model/Summary/Selectivity_vs_modelParams'];
    saveFig(gcf,outP,[{'tiff'} {'pdf'}])
end












