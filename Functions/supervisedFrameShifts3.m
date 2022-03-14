function shifts = supervisedFrameShifts3(meanFrames)
    
    numROIs = 1;
    shifts = zeros(3,length(meanFrames));
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
        imagesc(-ref) %(filter2(hSmall,meanFrames{mfi}) - filter2(hLarge, meanFrames{mfi}))
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
    
    minOverlap = nanmin(nansum(nansum(~isnan(cat(3,roiMF{:})),1),2))./2;
    xcLims = 60;
    rotLims = 10;
    allXC = nan(xcLims.*2+1,xcLims.*2+1,rotLims.*2+1,length(roiMF),length(roiMF));
    for i = round(length(roiMF)./2)
        for j = [1:i-1 i+1:length(roiMF)]
            allXC(:,:,:,i,j) = pvxcorr3(roiMF{i},roiMF{j},[xcLims xcLims rotLims],minOverlap);
            [x y z] = find(allXC(:,:,:,i,j) == nanmax(nanmax(nanmax(allXC(:,:,:,i,j)))));
            ts = [x-ceil(length(allXC(:,:,1,1,1))./2) y-ceil(length(allXC(:,:,1,1,1))./2) ...
                z-ceil(length(allXC(1,1,:,1,1))./2)];
            shifts(:,j) = ts;
        end
    end
end