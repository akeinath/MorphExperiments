function [nmX1 nmY1] = warpMaps(m1,mask)
    
    if length(mask(1,1,:))==1
        mask = repmat(double(mask),[1 1 length(m1(1,1,:))]);
        mask(mask==0) = nan;
    end
    m2 = mask;

    sm1 = ~isnan(m1(:,:,1));
    sm2 = ~isnan(m2(:,:,1));
    
    [x1 y1] = meshgrid(1:length(m1(1,:,1)),1:length(m1(:,1,1)));
    [x2 y2] = meshgrid(1:length(m2(1,:,1)),1:length(m2(:,1,1)));
    
    tmp1 = [sm1.*x1];
    tmp1(tmp1==0) = nan;
    sx1 = [nanmin(tmp1,[],2) nanmax(tmp1,[],2)];
    
    tmp1 = [sm1.*y1];
    tmp1(tmp1==0) = nan;
    sy1 = [nanmin(tmp1,[],1)' nanmax(tmp1,[],1)'];
    
    tmp1 = [sm2.*x2];
    tmp1(tmp1==0) = nan;
    sx2 = [nanmin(tmp1,[],2) nanmax(tmp1,[],2)];
    
    tmp1 = [sm2.*y2];
    tmp1(tmp1==0) = nan;
    sy2 = [nanmin(tmp1,[],1)' nanmax(tmp1,[],1)'];
    
    warpX = [nanmax(sx1(:,1),sx2(:,1)) nanmin(sx1(:,2),sx2(:,2))];
    warpY = [nanmax(sy1(:,1),sy2(:,1)) nanmin(sy1(:,2),sy2(:,2))];
    
    warpX = sort(warpX,2);
    warpY = sort(warpY,2);
    
    nmX1 = nan(size(m1));
    nmX2 = nan(size(m2));
    nmY1 = nan(size(m1));
    nmY2 = nan(size(m2));
    power = 2;
    for i = 1:length(warpX(:,1))
        if isnan(warpX(i,1))
            continue
        end
        
        try
            nmX1(i,warpX(i,1):warpX(i,2),:) = interp1(...
                [[([1:nansum(sm1(i,:))]-1)./nansum(sm1(i,:))].*(warpX(i,2)-warpX(i,1))]',...
                permute(m1(i,sm1(i,:),:),[3 2 1])',[0:warpX(i,2)-warpX(i,1)]','pchip','extrap');
        end
        try
            nmX2(i,warpX(i,1):warpX(i,2),:) = interp1(...
                [[([1:nansum(sm2(i,:))]-1)./nansum(sm2(i,:))].*(warpX(i,2)-warpX(i,1))]',...
                permute(m2(i,sm2(i,:),:),[3 2 1])',[0:warpX(i,2)-warpX(i,1)]','pchip','extrap');   
        end
    end
    
    for i = 1:length(warpY(:,1))
        if isnan(warpY(i,1))
            continue
        end
        try
            nmY1(warpY(i,1):warpY(i,2),i,:) = interp1(...
                [[([1:nansum(sm1(:,i))]-1)./nansum(sm1(:,i))].*(warpY(i,2)-warpY(i,1))]',...
                permute(m1(sm1(:,i),i,:),[3 1 2])',[0:warpY(i,2)-warpY(i,1)]','pchip','extrap');
        end
        try
            nmY2(warpY(i,1):warpY(i,2),i,:) = interp1(...
                [[([1:nansum(sm2(:,i))]-1)./nansum(sm2(:,i))].*(warpY(i,2)-warpY(i,1))]',...
                permute(m2(sm2(:,i),i,:),[3 1 2])',[0:warpY(i,2)-warpY(i,1)]','pchip','extrap');
        end
             
    end
    
    nmC1 = nanmean(cat(4,nmX1,nmY1),4);
    nmC2 = nanmean(cat(4,nmX2,nmY2),4);
end