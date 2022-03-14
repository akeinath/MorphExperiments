function [angdiff driftVariance driftAmount stressValue] = mds3D(mat,envs,envLabel,root)

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

    
%     [mdssim stress] = mdscale(tmp,3);
    
    [mdssim stressValue] = mdscale(tmp,3,'criterion','stress','start','cmdscale');
    
%     eigenVals./nanmax(eigenVals);
%     fprintf(['\n\tMax relative error (1D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (2D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:2))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (3D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:3))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (4D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:4))))) / max(tmp))]);
    mdssim = mdssim(:,1:3);
    
    figure
    set(gcf,'position',[50 50 800 400],'color','w')
    hold on
    keyShape = ['oooooo'];
    keyColor = envcm;
    initSize = 6;
    keySize = initSize+[floor([1:length(mat)]./6)].*2;
    lim = nanmax(abs(mdssim(:)))+nanmax(abs(mdssim(:))).*0.1;
    text(lim,lim,['Dim. 2'],'fontname','arial',...
        'fontsize',12,'horizontalalignment','center','verticalalignment','bottom')
    text(3.*lim,lim,['Dim. 3'],'fontname','arial',...
        'fontsize',12,'horizontalalignment','center','verticalalignment','bottom')
    text(lim.*2,[0],['Dim. 1'],'fontname','arial',...
        'fontsize',12,'horizontalalignment','center','verticalalignment','bottom','rotation',90)
    text(2.*lim,lim,['Stress: ' sprintf('%0.3f',stressValue)],'fontname','arial',...
        'fontsize',12,'horizontalalignment','center','verticalalignment','bottom')
    plot([0 lim.*4],[lim lim],'color','k')
    plot([lim.*2 lim.*2],[-lim lim],'color','k')
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        hold on
        plot((0).*2.*lim+lim+mdssim(envI,2),mdssim(envI,1),'marker','none', ...
            'color',keyColor(i,:),'linewidth',1);
        plot((1).*2.*lim+lim+mdssim(envI,3),mdssim(envI,1),'marker','none', ...
            'color',keyColor(i,:),'linewidth',1);
    end
    for i = 1:length(mdssim(:,1))
        envI = find(ismember(envLabel,envs(i)));
        h1 = plot((0).*2.*lim+lim+mdssim(i,2),mdssim(i,1),'marker',keyShape(envI), ...
            'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
            'markeredgecolor','w');
        h2 = plot((1).*2.*lim+lim+mdssim(i,3),mdssim(i,1),'marker',keyShape(envI), ...
            'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
            'markeredgecolor','w');
        if i == nanmax(find(ismember(envs,envLabel(envI))))
            set(h1,'markeredgecolor','k','linewidth',1.5);
            set(h2,'markeredgecolor','k','linewidth',1.5);
        end
    end
    axis off
    checkP(root);
    saveFig(gcf,[root],[{'pdf'} {'tiff'}]);

    
    figure
    set(gcf,'position',[50 50 400 400],'color','w')
    hold on
    lim = nanmax(abs(mdssim(:)))+0.1;
    set(gca,'xlim',[-lim lim],'ylim',[-lim lim],'zlim',[-lim lim],'linewidth',1.5)
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        plot3(mdssim(envI,1),mdssim(envI,2),mdssim(envI,3),'marker','none', ...
            'color',keyColor(i,:),'linewidth',1);
    end

    allH = [];
    for i = 1:length(mdssim(:,1))
        envI = find(ismember(envLabel,envs(i)));
        h = plot3(mdssim(i,1),mdssim(i,2),mdssim(i,3),'marker',keyShape(envI), ...
            'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
            'markeredgecolor','w');
        if i == nanmax(find(ismember(envs,envLabel(envI))))
            allH(envI) = h;
            set(h,'markeredgecolor','k','linewidth',1.5);
%                 text(mdssim(i,1),mdssim(i,2),upper(envLabel(envI)),'fontname','arial',...
%                     'fontsize',6,'fontweight','bold','horizontalalignment','center',...
%                     'verticalalignment','middle')
        end
    end
    grid on
    zlabel('MDS Dim 3','fontname','arial','fontsize',14);
    ylabel('MDS Dim 2','fontname','arial','fontsize',14);
    xlabel('MDS Dim 1','fontname','arial','fontsize',14);
    for i = setxor(0:30:150,[0 90 180])
        for j = setxor(0:30:150,[0 90 180])
            view([i j]);
            saveFig(gcf,[root '/Multiview_' num2str(i) '_' num2str(j)],[{'pdf'} {'tiff'}]);
        end
    end
    
%     zlabel('MDS Dim 3');
%     ylabel('MDS Dim 2');
%     xlabel('MDS Dim 1');
% %         axis equal
%     grid on
%     clear frames
%     for i = 1:1:360
%         view(i+0.5,20)
%         drawnow
%         frames(i) = getframe(gcf);
% 
%         im = frame2im(frames(i)); 
%         [imind,cm] = rgb2ind(im,256); 
%         if i == 1
%             checkP(root);
%             imwrite(imind,cm,root,'gif', 'Loopcount',inf,'DelayTime',0.025); 
%         else
%             imwrite(imind,cm,root,'gif','DelayTime',0.025,'WriteMode','append'); 
%         end
%     end
end