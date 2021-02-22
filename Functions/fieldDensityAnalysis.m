function fieldDensityAnalysis(allAvgMaps,root)
    nsims = 1000;
    figure
    set(gcf,'position',[50 50 900 250])
    adm = nan(15,15,3);
    sadm = nan(15,15,3,nsims);
    for i = 1:3
        adm(:,:,i) = help_norm(allAvgMaps{i}(1:15,1:15,:));
    end
    
    bounds = [0 length(allAvgMaps{1}(1,1,:)) ...
        length(allAvgMaps{1}(1,1,:))+length(allAvgMaps{2}(1,1,:)) ...
        length(allAvgMaps{1}(1,1,:))+length(allAvgMaps{2}(1,1,:))+length(allAvgMaps{3}(1,1,:))];
    
    pool = cat(3,allAvgMaps{:});
    for si = 1:nsims
        pool = pool(:,:,randperm(length(pool(1,1,:))));
        for i = 1:3
            sadm(:,:,i,si) = help_norm(pool(1:15,1:15,bounds(i)+1:bounds(i+1)));
        end
    end
    
    for i = 1:3
        upper_pval = 1-nanmean(bsxfun(@gt,adm(:,:,i),sadm(:,:,i,:)),4);
        lower_pval = 1-nanmean(bsxfun(@lt,adm(:,:,i),sadm(:,:,i,:)),4);
        sig_pixels = (upper_pval < 0.025) | (lower_pval < 0.025); 
        subplot(1,3,i)
        imagesc(adm(:,:,i))
        axis equal
        axis off
        caxis([0 0.2])
        
        [x y] = find(lower_pval < 0.025);
        text(y,x,'*','horizontalalignment','center','verticalalignment','middle',...
            'fontname','arial','fontsize',7,'color','w')
        [x y] = find(upper_pval < 0.025);
        text(y,x,'*','horizontalalignment','center','verticalalignment','middle',...
            'fontname','arial','fontsize',7,'color','k')
        colorbar
%         subplot(2,3,i+3)
%         imagesc([upper_pval < 0.025]-[lower_pval < 0.025])
%         alpha(double(sig_pixels))
%         axis equal
%         axis off
%         colormap hot
    end
    
    saveFig(gcf,[root '/Overall/AverageMapAnalysis'],[{'pdf'} {'tiff'}])
end

function tmp = help_norm(tmp)
    normer = repmat(nanmax(nanmax(tmp,[],1),[],2),size(tmp(:,:,1)));
    tmp = tmp./normer;
    tmp = tmp>0.75;
    
    tmp = nanmean(tmp,3);
end