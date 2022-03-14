function SFPs = compCropSFPs(usfps,cropSize)
    bsfps = nansum((usfps>0&~isnan(usfps)),4);
    bsfps = bsfps./repmat(nansum(nansum(bsfps,1),2),[size(bsfps(:,:,1)) 1]);
    [x y] = meshgrid(1:length(bsfps(1,:,1,1)),1:length(bsfps(:,1,1,1)));
    
    cx = nansum(nansum(bsxfun(@times,bsfps,x),1),2);
    cy = nansum(nansum(bsxfun(@times,bsfps,y),1),2);
    c = round([permute(cy,[3 1 2]) permute(cx,[3 1 2])]);
    
    clear bsfps
    ps = ceil(nanmax(cropSize)./2);
    SFPs = nan([cropSize length(usfps(1,1,1,:)) length(c(:,1))]);
    for k = 1:length(c(:,1))
        doX = c(k,1)-(ps-1):c(k,1)+(ps-1);
        doY = c(k,2)-(ps-1):c(k,2)+(ps-1);
        goodX = doX>0 & doX<length(usfps(:,1,1,1));
        goodY = doY>0 & doY<length(usfps(1,:,1,1));
        SFPs(goodX,goodY,:,k) = permute(usfps(doX(goodX),doY(goodY),k,:),[1 2 4 3]);
    end
end