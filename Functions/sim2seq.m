function [morphPoint maxDecorr] = sim2seq(sim,uIsPC,envs,doComps,root)
    if isempty(uIsPC)
        uIsPC = true([length(sim(1,1,:)) length(sim(:,1,1))]);
    end

    sigParams = [];
    morphPoint = [];
    maxDecorr = [];
    
    figure
    set(gcf,'position',[50 50 300.*floor(length(envs)./6) 250])
    for i = 1:floor(length(envs)./6)
        goodS = ((i-1).*6)+1:nanmin(((i-1).*6)+8,length(envs));
        tv = sim(goodS,goodS,:);
        tpc = uIsPC(:,goodS);
        ttv = repmat({[]},[size(tv(:,:,1))]);
        for j = 1:length(ttv)
            for k = 1:length(ttv)
                ttv{j,k} = permute(tv(j,k,tpc(:,j)|tpc(:,k)),[3 1 2]);
            end
        end
        subplot(1,floor(length(envs)./6),i)
        params = transitionPlot(ttv,envs(goodS),doComps,[]);
        if ~isempty(params)
            morphPoint = [morphPoint (params.sigmoidal_intercept)];
            maxDecorr = [maxDecorr params.sigmoidal_max_sep];
        end
    end
    saveFig(gcf,[root '/SequenceAnalyses_All'],[{'tiff'} {'pdf'}]);
    close all
    drawnow
end