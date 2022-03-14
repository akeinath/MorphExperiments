function [h doCells] = cumHist(d,x,doPlot)

    if nargin < 3 || isempty(doPlot)
        doPlot = true;
    end

    h = [];

    if isnumeric(d)
        tmp = [];
        for i = 1:length(d(1,:))
            tmp = [tmp {d(:,i)}];
        end
        d = tmp;
    end

    if doPlot
        set(gca,'xlim',[x(1) x(end)])
        hold on
        colors = [[0.1:0.1:1]' [0.1:0.1:1]' [1:-0.1:0.1]'];
    end
    
    abins = [];
    at = [];
    
    colors = circshift([[0.5:1./length(d)./2:1]' [0.5:1./length(d)./2:1]' [1:-1./length(d)./2:0.5]'],[0 -1]);
    ls = repmat({'-'},[1 length(d)]);
    for i = 1:length(d)
        if isempty(d{i})
            continue
        end
        [t bin] = histc(d{i},x);
        abins = [abins {bin}];
        t(end) = [];
%         t = cumsum(t);

        at = [at t];
% % %         t = t./(nansum(t));
%         h(i) = plot([x(1) x(1:end-1)+((x(2)-x(1))/2) x(end)],[0 t(:)' t(end)],'linewidth',2,...
%             'color',colors(i,:).*0.85,'linestyle',ls{i});
        t = [0; t];
        
        if doPlot
            plot([x(1:end-1); x(1:end-1); x(2:end)], ...
                [t(1:end-1)'; t(2:end)'; t(2:end)'],'linewidth',2,...
                'color',colors(i,:).*0.85,'linestyle',ls{i});
            hold on
        end
    end
    
    % create matched histogram indices
    npb = nanmin(at,[],2);
    doCells = repmat({[]},[1 length(d)]);
    for i = 1:length(npb)
        for j = 1:length(d)
            ci = find(abins{j}==i);
            ci = ci(randperm(length(ci)));
            doCells{j} = [doCells{j}; ci(1:npb(i))];
        end
    end
    doCells = cat(2,doCells{:});
end


































