function [allPrepped shifts] = supervisedFrameShifts(meanFrames,allPrepped)
    
    refSession = round(length(meanFrames)./2);

    numROIs = 1;
    roiMF = meanFrames;
    for mfi = 1:length(meanFrames)
        figure(1)
        set(gcf,'position',[50 50 800 800])
        
        ref = normMeanFrame(meanFrames{mfi},3);
        [x y] = meshgrid([1:length(ref(1,:,1))]-length(ref(1,:,1))./2+0.5, ...
        [1:length(ref(:,1,1))]-length(ref(:,1,1))./2+0.5);
        d2c = sqrt(x.^2 + y.^2);
        mask = d2c < nanmin(size(ref(:,:,1))).*0.5;
        ref(~mask) = nan;
        imagesc(ref) %(filter2(hSmall,meanFrames{mfi}) - filter2(hLarge, meanFrames{mfi}))
        colormap gray
        axis equal
        
        mask = false(size(ref));
        for i = 1:numROIs
            rect = getrect(); 
            roi = uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
            mask(roi(3):roi(4),roi(1):roi(2)) = true;
        end
        
        ref(~mask) = nan;
        roiMF{mfi} = ref;
    end
    close all
    drawnow
    
    scaleFactors = [1];
    
    fprintf('\n\t\t\tAligning mean frames by crosscorrelation...  ');
    tic
    xcLims = 70;
    rotLim = 10;
    shifts = nan(3,length(roiMF),length(scaleFactors));
    for sfi = 1:length(scaleFactors)
        rroiMF = imresize(cat(3,roiMF{:}),scaleFactors(sfi));
        if sfi > 1
            
            preshift = nansum(bsxfun(@times,scaleFactors(sfi)./repmat(permute(scaleFactors(1:sfi-1),[1 3 2]),[sfi-1 1 1]), ...
                shifts(1:2,:,1:sfi-1)),3);
            
            for i = 1:length(rroiMF(1,1,:))
                tmp = rroiMF(:,:,i);
                ps = xcLims;
                tmp = padarray(tmp,[ps ps],nan);
                tmp = circshift(tmp,[preshift(1:2,i)']);
                tmp = tmp(ps+1:end-ps,ps+1:end-ps,:);
                rroiMF(:,:,i) = tmp;
            end
        end
        minOverlap = nanmin(nansum(nansum(~isnan(rroiMF),1),2))./2;
        for i = refSession
            for j = [1:i-1 i+1:length(roiMF)]
                if sfi == length(scaleFactors)
                    tmp = pxcorr(rroiMF(:,:,i),rroiMF(:,:,j),[xcLims xcLims rotLim],minOverlap);
                else
                    tmp = pxcorr(rroiMF(:,:,i),rroiMF(:,:,j),[xcLims xcLims 0],minOverlap);
                end
                x = find(tmp == nanmax(nanmax(nanmax(tmp))));
                [x y z] = ind2sub(size(tmp(:,:,:)),x);
                ts = [nanmedian(x)-ceil(length(tmp(:,:,1))./2) ...
                    nanmedian(y)-ceil(length(tmp(:,:,1))./2) ...
                    nanmedian(z)-ceil(length(tmp(1,1,:))./2)];
                shifts(:,j,sfi) = ts;
            end
        end
    end
    totalShift = nansum(bsxfun(@times,scaleFactors(end)./repmat(permute(scaleFactors,[1 3 2]),[sfi-1 1 1]), ...
        shifts(1:2,:,1:sfi)),3); 
    
    ps = xcLims./nanmin(scaleFactors);
    for i = setdiff(1:length(allPrepped),refSession)
        tmp = allPrepped{i};
        tmp = permute(tmp,[2 3 1]);
        parfor k = 1:length(tmp)
            tmp2 = imresize(tmp(:,:,k),scaleFactors(end));
            tmp2 = padarray(tmp2,[ps ps],0);
            tmp2 = circshift(tmp2,[totalShift(1:2,i)']);
            tmp2 = tmp2(ps+1:end-ps,ps+1:end-ps,:);
            tmp2 = imresize(tmp2,1./scaleFactors(end));
            tmp(:,:,k) = tmp2;
        end
        allPrepped{i} = permute(tmp,[3 1 2]);
    end
    
    tmp = toc;
    fprintf('  %0.3fs.\n',tmp);
    close all
    drawnow;
end



















