function [actual actual_sim] = contextSVM(maps,envs,nDays)
    close all
    drawnow

    trainSets = bsxfun(@plus,1:nDays,[1:32-nDays+1]'-1);    
    
    [actual actual_sim] = help_contextClassification(maps,envs,trainSets(1:6:end,:),false,false);
    
%     [cross_actual cross_actual_sim] = help_contextClassification(maps,envs,trainSets(1:1:end,:),false,true);
    
%     cross_actual = help_contextClassification(maps,envs,trainSets,false,true);
%     
%     nsims = 20;
%     null = nan(nsims,length(actual));
%     tic
%     for si = 1:20
%         null(si,:) = help_contextClassification(maps,envs,trainSets,true);
%     end
%     toc
    
%     close all
%     drawnow
%     figure
%     set(gcf,'position',[50 50 400 250])
%     plot(1:length(actual(1,:)),actual,'color',[0.9 0.6 0.3])
%     hold on
%     plot(get(gca,'xlim'),[0.5 0.5],'color',[0.5 0.5 0.5],'linestyle','--')
%     set(gca,'ylim',[0 1])
%     ylabel('Prediction Accuracy')
 
end