function plotPaths(paths)
    clc
    close all
    drawnow
    

    warning off all
     %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    labels = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    
    pause_thresh = -2;
    envSize = [17 17];
    
    allPropAligned = [];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);

        
        figure
        set(gcf,'position',[50 50 8.*100 ceil(length(sessions)./8).*100])
        
        fprintf(['\t\tLoading Data... '])
        for si = 1:length(sessions)
            s = load(sessions{si},'processed');
            
            subplot(ceil(length(sessions)./8),8,si)
            plot(s.processed.p(2,1:5:end),s.processed.p(1,1:5:end),'color','k')
            axis square
            set(gca,'xlim',[-2 40],'ylim',[-2 40])
            axis off
        end  
        
        mouseName = upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end);
        root = ['Plots/Paths_Grouped/' mouseName];
        saveFig(gcf,root,[{'tiff'} {'pdf'} {'jpeg'}])
    end
end