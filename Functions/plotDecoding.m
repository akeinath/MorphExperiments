function plotDecoding(val,root)

    set(0,'DefaultAxesColorOrder',plasma(7))';

    figure
    set(gcf,'position',[50 50 750 225])
    subplot(1,3,1)
    plot(permute(val(:,2,:),[3 1 2])','marker','o','markerfacecolor','auto')
    set(gca,'xlim',[0 length(val(:,1,1))+1],'xtick',[1:length(val(:,1,1))], ...
        'xticklabel',[{'1-6'} {'7-12'} {'13-18'} {'19+'}])
    hold on
    plot(nanmean(permute(val(:,2,:),[3 1 2]),1)','marker','o','markerfacecolor','k','color','k')
    plot(get(gca,'xlim'),[(1./2) (1./2)],'color','k','linestyle','--')
    set(gca,'ylim',[0 1])
    subplot(1,3,2)
    plot(permute(val(:,1,:),[3 1 2])','marker','o','markerfacecolor','auto')
    set(gca,'xlim',[0 length(val(:,1,1))+1],'xtick',[1:length(val(:,1,1))], ...
        'xticklabel',[{'1-6'} {'7-12'} {'13-18'} {'19+'}])
    hold on
    plot(nanmean(permute(val(:,1,:),[3 1 2]),1)','marker','o','markerfacecolor','k','color','k')
    plot(get(gca,'xlim'),[(1./6) (1./6)],'color','k','linestyle','--')
    set(gca,'ylim',[0 0.6])
    subplot(1,3,3)
    plot(permute(val(:,3,:),[3 1 2])','marker','o','markerfacecolor','auto')
    set(gca,'xlim',[0 length(val(:,1,1))+1],'xtick',[1:length(val(:,1,1))], ...
        'xticklabel',[{'1-6'} {'7-12'} {'13-18'} {'19+'}])
    hold on
    plot(nanmean(permute(val(:,3,:),[3 1 2]),1)','marker','o','markerfacecolor','k','color','k')
    plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
    set(gca,'ylim',[-0.2 1])
    saveFig(gcf,root,[{'tiff'} {'pdf'}])
    
    fid = fopen([root '_stats.txt'],'w');
    fprintf(fid,'\n\n\t\t\t%%%%%%%%%% Context Group Prediction %%%%%%%%%%\n');
    [a b c d] = ttest(permute(val(:,2,:),[3 1 2]),(1./2),'tail','right');
    fprintf(fid,'\n\t t(%i) = %0.2f, p = %.2e',[d.df; d.tstat; b]);
    [pval,tbl] = anova1(permute(val(:,2,:),[3 1 2]),[],'off');
    fprintf(fid,'\n\t F(%i,%i) = %0.2f, p = %.2e',[tbl{2,3}; tbl{3,3}; tbl{2,5}; tbl{2,6}]);
    fprintf(fid,'\n\n\t\t\t%%%%%%%%%% Environment ID Prediction %%%%%%%%%%\n');
    [a b c d] = ttest(permute(val(:,1,:),[3 1 2]),(1./6),'tail','right');
    fprintf(fid,'\n\t t(%i) = %0.2f, p = %.2e',[d.df; d.tstat; b]);
    [pval,tbl] = anova1(permute(val(:,1,:),[3 1 2]),[],'off');
    fprintf(fid,'\n\t F(%i,%i) = %0.2f, p = %.2e',[tbl{2,3}; tbl{3,3}; tbl{2,5}; tbl{2,6}]);
    fprintf(fid,'\n\n\t\t\t%%%%%%%%%% Best Match Similarity %%%%%%%%%%\n');
    [pval,tbl] = anova1(permute(val(:,3,:),[3 1 2]),[],'off');
    fprintf(fid,'\n\t F(%i,%i) = %0.2f, p = %.2e',[tbl{2,3}; tbl{3,3}; tbl{2,5}; tbl{2,6}]);
    tmp = permute(val(:,3,:),[3 1 2]);
    for i = 1:length(tmp(1,:))
        for j = i+1:length(tmp(1,:))
            [a b c d] = ttest(tmp(:,i),tmp(:,j));
            fprintf(fid,'\n\t Comparison(%i, %i) t(%i) = %0.2f, p = %.2e',[i; j; d.df; d.tstat; b]);
        end
    end
    fclose(fid);
    
%     xlswrite([root '_ContextGroupVals.xls'],permute(val(:,2,:),[3 1 2]));
%     xlswrite([root '_EnvironmentIDVals.xls'],permute(val(:,2,:),[3 1 2]));
end




























