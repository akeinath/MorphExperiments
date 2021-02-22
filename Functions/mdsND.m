function mdsND(mat,envs,envLabel,root,nDims)

    if nargin < 5 || isempty(nDims)
        nDims = 3;
    end

    %%% Multidimensional 3D video
%     tmp = cellfun(@nanmedian,kCrossSim);
    if iscell(mat)
        tmp = cellfun(@nanmedian,mat);
    else
        tmp = nanmedian(mat,3);
    end
    tmp = nanmax(tmp,tmp');
    tmp(logical(eye(size(tmp)))) = 1;
    tmp = 1-tmp;

    
%     [mdssim stress] = mdscale(tmp,nDims);
    
    [mdssim eigenVals] = cmdscale(tmp);
    
%     eigenVals./nanmax(eigenVals);
%     fprintf(['\n\tMax relative error (1D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (2D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:2))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (3D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:3))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (4D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:4))))) / max(tmp))]);
    mdssim = mdssim(:,1:nDims);

    for combI = 1:nDims
        for combJ = combI+1:nDims
    
            q = figure;
            set(gcf,'position',[50 50 400 400])
            lim = nanmax(abs(mdssim(:)))+0.1;
            set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
            hold on
            keyShape = ['sssooo'];
            keyColor = [0.9 0.1 0.1; 0.9 0.4 0.4; 0.9 0.7 0.7; 0.7 0.7 0.9; 0.4 0.4 0.9; 0.1 0.1 0.9].*0.5+0.5;
            initSize = 6;
            keySize = initSize+[floor([1:length(tmp)]./6)].*2.5;
            for i = 1:length(envLabel)
                envI = find(ismember(envs,envLabel(i)));
                plot(mdssim(envI,combI),mdssim(envI,combJ),'marker','none', ...
                    'color',keyColor(i,:),'linewidth',1);
            end

            for i = 1:length(mdssim(:,1))
                envI = find(ismember(envLabel,envs(i)));
                h = plot(mdssim(i,combI),mdssim(i,combJ),'marker',keyShape(envI), ...
                    'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
                    'markeredgecolor','w');
                if i == nanmax(find(ismember(envs,envLabel(envI))))
                    set(h,'markeredgecolor','k','linewidth',1.5);
                    text(mdssim(i,combI),mdssim(i,combJ),upper(envLabel(envI)),'fontname','arial',...
                        'fontsize',6,'fontweight','bold','horizontalalignment','center',...
                        'verticalalignment','middle')
                end
            end
            ylabel(['MDS Dim ' num2str(combJ)]);
            xlabel(['MDS Dim ' num2str(combI)]);
            axis square
            
            
            checkP([root 'Projection' num2str(combI) '-' num2str(combJ)]);
            saveFig(gcf,[root 'Projection' num2str(combI) '-' num2str(combJ)],[{'pdf'} {'tiff'}]);
            close(q);
            drawnow;
        end
    end
end