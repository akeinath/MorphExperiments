function pred = DO2S(x,p)
    
    [mx my] = meshgrid(x,x);
    mx = (mx-nanmean(mx(:)));
    my = (my-nanmean(my(:)));
    
    rotMat1 = [cosd(45) sind(45); -sind(45) cosd(45)];
    t1 = [1./(1 + exp(-p(2).*([1 0]*rotMat1*[mx(:) my(:)]') + p(1)))].*p(3) + p(4);
    pred1 = reshape(t1,size(mx));
    
    rotMat2 = [cosd(45) sind(45); -sind(45) cosd(45)];
    t2 = [1./(1 + exp(-p(6).* abs([0 1]*rotMat2*[mx(:) my(:)]') + p(5)))].*p(7) + p(8);
    pred2 = reshape(t2,size(mx));
    
    pred = bsxfun(@times,pred1,pred2'); %+p(end);
end