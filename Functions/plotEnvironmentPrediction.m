function [binned_out binned_sim_out]= plotEnvironmentPrediction(val,sim,root,maxDays)
    
    if nargin < 3 || isempty(maxDays)
        maxDays = length(across.contextPrediction.actual(1,:,1));
    end
    
    toPlot = [];
    for i = 1:length(val(1,1,:))
        toPlot = [toPlot {val(:,1:maxDays,i)}];
    end
    
    doC = (1-v2rgb(1:length(toPlot)+1)).*0.5 + 0.25;
    doC(1,:) = [];
    
    figure
    set(gcf,'position',[50 50 400 300])
    if length(toPlot) == 1
        for i = 1:length(toPlot)
            plot(toPlot{i}','color',doC(i,:),'linestyle','-',...
                'linewidth',0.5);
        end
    end
    hold on
    h = mkLine(toPlot,[],doC);
    plot(get(gca,'xlim'),[1./6 1./6],'color',[0.5 0.5 0.5],'linestyle','--')
    set(gca,'ylim',[0 1])
    ylabel('Prediction Accuracy')
    xlabel('Lag from training set')
    legend(h(1,:),[{'1st'} {'2nd'} {'3rd'} {'4th'}],'location','southeast')
    saveFig(gcf,[root 'EnvironmentPrediction_Accuracy'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    batchSize = 6;
    toPlot = [];
    for i = 1:length(val(1,1,:))
        for j = 1:maxDays./batchSize
            if j == nanmax(maxDays./batchSize)
                toPlot{j,i} = [nanmean(val(:,(j-1).*batchSize+1:end,i),2)];
            else
                toPlot{j,i} = [nanmean(val(:,(j-1).*batchSize+1:(j).*batchSize,i),2)];
            end
        end
    end
    
    figure
    set(gcf,'position',[50 50 250 225])
    mkBar(toPlot',[{'1-6'} {'7-12'} {'13-18'} {'19+'}],doC)
    hold on
    plot(get(gca,'xlim'),[1./6 1./6],'color','k','linestyle','--')
    saveFig(gcf,[root 'EnvironmentPrediction_Accuracy_Bar'],[{'tiff'} {'pdf'} {'jpeg'}])
    binned_out = toPlot;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    toPlot = [];
    for i = 1:length(sim(1,1,:))
        toPlot = [toPlot {sim(:,1:maxDays,i)}];
    end
    
    figure
    set(gcf,'position',[50 50 400 300])
    if length(toPlot) == 1
        for i = 1:length(toPlot)
            plot(toPlot{i}','color',doC(i,:),'linestyle','-',...
                'linewidth',0.5);
        end
    end
    hold on
    h = mkLine(toPlot,[],doC);
    plot(get(gca,'xlim'),[0.0 0.0],'color',[0.5 0.5 0.5],'linestyle','--')
    set(gca,'ylim',[-0.2 1])
    ylabel('Max PV Correlation (r)')
    xlabel('Lag from training set')
    legend(h(1,:),[{'1st'} {'2nd'} {'3rd'} {'4th'}],'location','southeast')
    saveFig(gcf,[root 'EnvironmentPrediction_Similarity'],[{'tiff'} {'pdf'} {'jpeg'}])
    
    
    
    
    batchSize = 6;
    toPlot = [];
    for i = 1:length(sim(1,1,:))
        for j = 1:maxDays./batchSize
            if j == nanmax(maxDays./batchSize)
                toPlot{j,i} = [nanmean(sim(:,(j-1).*batchSize+1:end,i),2)];
            else
                toPlot{j,i} = [nanmean(sim(:,(j-1).*batchSize+1:(j).*batchSize,i),2)];
            end
        end
    end
    binned_sim_out = toPlot;
end






























