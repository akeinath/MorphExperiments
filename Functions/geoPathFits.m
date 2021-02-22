function vals = geoPathFits(p1,p2,gT1,gT2)
    vals = [];
    
    m(:,:,:,1) = mkTraceMaps(p1,gT1,[],[17 17]);
    m(:,:,:,5) = mkTraceMaps(p2,gT2,[],[17 17]);
    
    
    [p1x p2x p1y p2y] = fitDeformationPaths(p1,p2);
    
    m(:,:,:,2) = mkTraceMaps(p1x,gT1,[],[17 17]);
    m(:,:,:,6) = mkTraceMaps(p2x,gT2,[],[17 17]);
    m(:,:,:,3) = mkTraceMaps(p1y,gT1,[],[17 17]);
    m(:,:,:,7) = mkTraceMaps(p2y,gT2,[],[17 17]);
    m(:,:,:,4) = nanmean(m(:,:,:,[2 3]),4);
    m(:,:,:,8) = nanmean(m(:,:,:,[6 7]),4);
    
    for k = 1:length(m(1,1,:,1))
        
        figure
        set(gcf,'position',[50 50 800 400])
        for j = 1:8
            subplot(2,4,mod(j+3,8)+1)
            imagesc(m(:,:,k,j))
            alpha(double(~isnan(m(:,:,k,j))))
            axis off
            axis equal
            colormap jet
        end
        saveFig(gcf,['Plots/GeoFitExamples/Cell_' num2str(k)],[{'tiff'} {'pdf'}])
        close all
    end
    
    vals = nan([4 4 length(m(1,1,:,1))]);
    for i = 1:4
        for j = 5:8
            tm1 = m(:,:,:,i);
            tm2 = m(:,:,:,j);
            rm1 = reshape(tm1,[prod(size(m(:,:,1,1))) length(m(1,1,:,1))]);
            rm2 = reshape(tm2,[prod(size(m(:,:,1,1))) length(m(1,1,:,1))]);
            isBad = any(isnan(rm1),2) | any(isnan(rm2),2);
            tv = corr(rm1(~isBad,:),rm2(~isBad,:));
            vals(i,j-4,:) = tv(logical(eye(size(tv))));
        end
    end
    vals = permute(nanmax(nanmax(vals,[],1),[],2),[3 1 2]);
end