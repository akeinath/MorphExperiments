function [angdiff driftVariance driftAmount stressValue] = mds2D(mat,envs,envLabel,root,doPlot)
    
    if nargin < 5 || isempty(doPlot)
        doPlot = true;
    end

    %%% Multidimensional Scaling 2D
    if iscell(mat)
        tmp = cellfun(@nanmean,mat);
    else
        tmp = nanmean(mat,3);
    end
    tmp = nanmax(tmp,tmp');
    
    if doPlot
        figure()
        imagesc(tmp)
        caxis([0 0.8])
        alpha(double(~isnan(tmp)));
        colormap parula
        axis square
        axis off
        if nargin > 3 && ~isempty(root)
            checkP(root);
            saveFig(gcf,[root '_RDM'],[{'pdf'} {'tiff'}]);
        end
    end
    
    tmp(logical(eye(size(tmp)))) = 1;
    dm = 1-(tmp);
%     dm = sqrt(2.*(1-tmp));

    [mdssim stressValue] = mdscale(dm,2,'criterion','stress','start','cmdscale');

%     eigenVals./nanmax(eigenVals);
%     fprintf(['\n\tMax relative error (1D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (2D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:2))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (3D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:3))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (4D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:4))))) / max(tmp))]);
    mdssim = mdssim(:,1:2);
%     
%         factoran(mdssim,2);
    
    for i = 1:length(tmp)
        for j = i+1:length(tmp)
            
        end
    end


    features = [];
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        for k = 1:length(envI)
            features(envI(k),:) = [k i];
        end
    end
    
%     driftSlopes = [];
%     for i = 1:length(envLabel)
%         envI = find(ismember(envs,envLabel(i)));
%         driftSlopes = [driftSlopes; polyfit(mdssim(envI,1),mdssim(envI,2),1)];
%     end
%     
%     contextSlopes = [];
%     doS = 6;
%     for i = 1:doS:length(envs)-(doS)+1
%         contextSlopes = [contextSlopes; polyfit(mdssim(i:i+doS-1,1),mdssim(i:i+doS-1,2),1)];
%     end

    driftSlopes = [];
    for i = [1 length(envLabel)]
        envI = find(ismember(envs,envLabel(i)));
        driftSlopes = [driftSlopes; diff(mdssim(envI([1 end]),1)) diff(mdssim(envI([1 end]),2))];
    end
%     driftSlopes = driftSlopes([1 end],:);
    [t d] = cart2pol(driftSlopes(:,2),driftSlopes(:,1));
    driftVariance = rad2deg(diff(t));
    allDriftSlopes = t;
    driftSlopes = circ_mean(t);
    
    
    ads = [];
    for i = [1:length(envLabel)]
        envI = find(ismember(envs,envLabel(i)));
        ads = [ads; [1:length(envI)-1]' diff(mdssim(envI([1:end]),1)) diff(mdssim(envI([1:end]),2))];
    end
    [t d] = cart2pol(ads(:,2),ads(:,3));
    driftAmount = [];
    for i = 1:nanmax(ads(:,1))
        driftAmount = [driftAmount; nanmean(d(ads(:,1)==i))];
    end
    
    
%     contextSlopes = [];
%     a1 = find(ismember(envs,{'sq1'}));
%     a2 = find(ismember(envs,{'g1'}));
%     for i = 1:nanmin(length(a1),length(a2))
%         contextSlopes = [contextSlopes; (mdssim(a1(i),2)-mdssim(a2(i),2)) ./ ...
%             (mdssim(a1(i),1)-mdssim(a2(i),1))];
%     end

    contextSlopes = [];
    a1 = find(ismember(envs,{'sq1'}));
    a2 = find(ismember(envs,{'g1'}));
    for i = 1:nanmin(length(a1),length(a2))
        contextSlopes = [contextSlopes; (mdssim(a1(i),1)-mdssim(a2(i),1)) ...
            (mdssim(a1(i),2)-mdssim(a2(i),2))];
    end
    [t d] = cart2pol(contextSlopes(:,2),contextSlopes(:,1));
    contextSlopes = circ_mean(t);

%     angdiff = rad2deg(diff(atan([nanmedian(driftSlopes(:,1)) nanmedian(contextSlopes(:,1))])));
    
    angdiff = rad2deg(driftSlopes-contextSlopes);
    
    if doPlot
    
        figure
        set(gcf,'position',[50 50 400 400])
    %     subplot(1,2,1)
%         lim = nanmax(abs(mdssim(:)))+0.1;
        lim = 0.7;
        set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
        hold on
    %     for i = -1:0.2:1
    %         plot(-1:0.01:1,polyval([driftPoly(:,1) i],-1:0.01:1),'color',[0.85 0.85 0.85], ...
    %             'linewidth',0.5,'linestyle','--')
    %     end
        keyShape = ['sssooo'];
    %     keyColor = [0.9 0.1 0.1; 0.9 0.4 0.4; 0.9 0.7 0.7; 0.7 0.7 0.9; 0.4 0.4 0.9; 0.1 0.1 0.9].*0.5+0.5;
        keyColor = envcm;
        initSize = 6;
        keySize = initSize+[floor([1:length(tmp)]./6)].*2.5;
        for i = 1:length(envLabel)
            envI = find(ismember(envs,envLabel(i)));
            plot(mdssim(envI,1),mdssim(envI,2),'marker','none', ...
                'color',keyColor(i,:),'linewidth',1);
        end

        for i = 1:length(mdssim(:,1))
            envI = find(ismember(envLabel,envs(i)));
            h = plot(mdssim(i,1),mdssim(i,2),'marker',keyShape(envI), ...
                'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
                'markeredgecolor','w');
            if i == nanmax(find(ismember(envs,envLabel(envI))))
                set(h,'markeredgecolor','k','linewidth',1.5);
    %             text(mdssim(i,1),mdssim(i,2),upper(envLabel(envI)),'fontname','arial',...
    %                 'fontsize',6,'fontweight','bold','horizontalalignment','center',...
    %                 'verticalalignment','middle')
            end
        end
        ylabel('MDS Dim 2 (a.u.)');
        xlabel('MDS Dim 1 (a.u.)');
        axis square

        hold on
%         as = atan([nanmedian(driftSlopes(:,1)) nanmedian(contextSlopes(:,1))]);
        h1 = plot([-1 1]*sin(driftSlopes(1)),[-1 1]*cos(driftSlopes(1)),'color','k','linestyle','-','linewidth',1);
        h2 = plot([-1 1]*sin(contextSlopes(1)),[-1 1]*cos(contextSlopes(1)),'color','k','linestyle','--','linewidth',1);

        %%% subtract drift
    %     for i = 1:length(envLabel)
    %         envI = find(ismember(envs,envLabel(i)));
    %         subsim(envI,2) = subsim(envI,2) - polyval(allSlopes(i,:),subsim(envI,1))
    %     end
    %     
    %     theta = cart2pol(1,driftPoly(:,1));
    %     rotmat = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
    %     subsim = mdssim*rotmat;

    %     subplot(1,2,2)
    %     lim = nanmax(abs(subsim(:)))+0.1;
    %     set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
    %     hold on
    %     keyShape = ['sssooo'];
    %     keyColor = [0.9 0.1 0.1; 0.9 0.4 0.4; 0.9 0.7 0.7; 0.7 0.7 0.9; 0.4 0.4 0.9; 0.1 0.1 0.9].*0.5+0.5;
    %     initSize = 6;
    %     keySize = initSize+[floor([1:length(tmp)]./6)].*2.5;
    %     for i = 1:length(envLabel)
    %         envI = find(ismember(envs,envLabel(i)));
    %         plot(subsim(envI,1),subsim(envI,2),'marker','none', ...
    %             'color',keyColor(i,:),'linewidth',1);
    %     end
    % 
    %     for i = 1:length(subsim(:,1))
    %         envI = find(ismember(envLabel,envs(i)));
    %         h = plot(subsim(i,1),subsim(i,2),'marker',keyShape(envI), ...
    %             'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
    %             'markeredgecolor','w');
    %         if i == nanmax(find(ismember(envs,envLabel(envI))))
    %             set(h,'markeredgecolor','k','linewidth',1.5);
    %             text(subsim(i,1),subsim(i,2),upper(envLabel(envI)),'fontname','arial',...
    %                 'fontsize',6,'fontweight','bold','horizontalalignment','center',...
    %                 'verticalalignment','middle')
    %         end
    %     end
    %     ylabel('MDS Dim 2');
    %     xlabel('MDS Dim 1');
    %     axis square

        text(0,lim,['Stress: ' sprintf('%0.3f',stressValue)],'fontname','arial',...
            'fontsize',12,'horizontalalignment','center','verticalalignment','bottom')
        if nargin > 3 && ~isempty(root)
            checkP(root);
            saveFig(gcf,[root],[{'pdf'} {'tiff'}]);
        end
    end
end