function pred = DOL(x,p)

    [mx my] = meshgrid(x,x);
    
%     rotMat1 = [cosd(45) sind(45); -sind(45) cosd(45)];
    rotMat1 = [cosd(0) sind(0); -sind(0) cosd(0)];
    rm = [mx(:) my(:)]*rotMat1;
    
    lin1 = polyval([p(1:(length(p)-1)./2) 0],rm(:,1));
    lin2 = polyval([p((length(p)-1)./2 +1:end-1) 0],rm(:,2));
    
    pred1 = reshape(lin1,size(mx));
    pred2 = reshape(lin2,size(mx));
    
    pred = pred1.*pred2+p(end);
    
%     pred = bsxfun(@minus,lin1,lin2')+p(end);
    
    mask = triu(true(size(pred)),1);
    tv = pred(mask);
    pred = pred';
    pred(mask) = tv;
end
