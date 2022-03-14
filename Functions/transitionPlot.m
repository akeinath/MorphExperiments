function params = transitionPlot(data,envs,doComps,root)
    params = [];

    if ~iscell(data)
        data = num2cell(data);
    end
    data(logical(eye(size(data)))) = {nan};
    
%     data = cellfun(@atanh,data,'uniformoutput',false);

    toPlot = repmat({[]},[1 length(doComps(:,1))]);
    compLabel = [];
    for i = 1:length(doComps(:,1))
        compLabel = [compLabel {[doComps{i,1} '-' doComps{i,2}]}];
        ga = bsxfun(@times,ismember(envs,doComps(i,1)),ismember(envs,doComps(i,2))');
        gb = bsxfun(@times,ismember(envs,doComps(i,2)),ismember(envs,doComps(i,1))');
        isGood = (ga|gb);
        if ~any(isGood)
            continue
        end
        inc = cat(1,data{isGood});
        toPlot{i} = inc(~isnan(inc));
    end
    
%     figure
%     set(gcf,'position',[50 50 25.*length(doComps(:,1)) 350])
    cm = envcm();
    mkAcrossLine([toPlot(1:2:end); toPlot(2:2:end)],[],cm([1 6],:));
%     mkGraph([toPlot(1:2:end); toPlot(2:2:end)],doComps(1:2:end,2)');
%     mkBow([toPlot(1:2:end); toPlot(2:2:end)]);
    hold on
    set(gca,'ylim',[0 0.8])
%     plot(get(gca,'xlim'),[0 0],'linewidth',1,'linestyle','--','color','k')
    xlabel('Condition Comparison')

%     figure
%     set(gcf,'position',[50 50 30.*length(doComps(:,1)) 350])
%     mkBow([toPlot(1:2:end); toPlot(2:2:end)],doComps(1:2:end,2)');
%     set(gca,'ylim',[-1 1])
%     hold on
%     plot(get(gca,'xlim'),[0 0],'linewidth',1,'linestyle','--','color','k')
%     xlabel('Condition Comparison')
    
    
%     %%% Analysis of curve
%     
    
    [sfit sint lfit mad max_sep] = fitSigLin([toPlot(1:2:end); toPlot(2:2:end)]);
    params = struct;
    params.sigmoidal_mad = mad;
    params.sigmoidal_max_sep = max_sep;
    params.sigmoidal_fit = sfit(:,1);
    params.sigmoidal_params = sfit(:,2:end);
    params.sigmoidal_intercept = sint;
    params.linear_fit = lfit;
    a = cellfun(@nanmean,toPlot(1:2:end));
    b = cellfun(@nanmean,toPlot(2:2:end));
    xa = 1:6;
    params.poly_intercept = polyxpoly(xa(~isnan(a)),a(~isnan(a)),xa(~isnan(b)),b(~isnan(b)));
    if length(params.poly_intercept)>1
        [qa qb] = nanmin(abs(params.poly_intercept-3.5));
        params.poly_intercept = params.poly_intercept(qb);
    end
    
%     plot(ones(1,2).*params.poly_intercept,get(gca,'ylim'),'color','k','linestyle','-.')
    if ~isempty(sint)
%         plot(sint,nanmax(get(gca,'ylim'))-0.1.*range(get(gca,'ylim')),'color','k', ...
%             'linestyle','none','marker','x','markersize',8,'linewidth',2)
    else
        params.poly_intercept = 3.5;
    end
    
    tmp = cellfun(@nanmean,toPlot);
    tmp = [tmp(1:2:end); tmp(2:2:end)];
    params.attractor_strength = nanmax(abs(tmp(1,:)-tmp(2,:)));
%     try
%         tmp([1 2],xa > params.poly_intercept) = tmp([2 1],xa > params.poly_intercept);
%         params.attractor_strength = nanmean(tmp(1,:)-tmp(2,:));
%     end
    
    %%% Save Figure
    if nargin > 3 && ~isempty(root)
        saveFig(gcf,root,[{'pdf'} {'tiff'}]);
    end
end