function interpNormCorreVid(ms)
%     clc
    
    fprintf('\n\t--Checking and correcting error frames...\t')
    tic
    
    obj = VideoReader([ms.dirName '/' 'msvideo.avi']);
    tmp = read(obj);
    tmp = double(tmp);
    
    doPix = true(numel(tmp(:,:,1,1)));
    pot = false(size(tmp(:,:,1,1)));
    [x y] = meshgrid(1:length(pot(1,:)),1:length(pot(:,1)));
    d2c = sqrt([x-round(length(x(1,:))./2)].^2 + [y-round(length(y(:,1))./2)].^2);
    isCenter = d2c < nanmean(size(pot)./4);
    doPix = isCenter(:);
    
    tmp2 = reshape(tmp,[obj.Width.*obj.Height length(tmp(1,1,1,:))]);
%     dmAct = diff(nanmean(tmp2));
    dmAct = nanmean(diff(tmp2(doPix,:),[],2));
        
    totalInterp = 0;
    thresh = nanstd(dmAct).*6; % toss 6 stdev above mean.
    
%     figure
%     set(gcf,'position',[50 50 600 300]);
%     plot(dmAct)
    while any(abs(dmAct) > thresh)        
%         plot(dmAct)
%         drawnow

        totalInterp = totalInterp+1;
        
        badInd = find(dmAct<-thresh);
        blah = diff(badInd);
        ignore = [false blah<600];
        badInd(ignore) = [];
        tmp(:,:,badInd+1) = tmp(:,:,badInd);
        
        badInd = find(dmAct>thresh);
        blah = diff(badInd);
        ignore = [blah<600 false];
        badInd(ignore) = [];
        tmp(:,:,badInd) = tmp(:,:,badInd-1);

        tmp2 = reshape(tmp,[obj.Width.*obj.Height length(tmp(1,1,1,:))]);
%         dmAct = diff(nanmean(tmp2));
        dmAct = nanmean(diff(tmp2(doPix,:),[],2));
        thresh = nanstd(dmAct).*6; % toss 4 stdev above mean.
    end
    fprintf(['\n\t\t--Number of interpolations: ' num2str(totalInterp) '\t'])
%     plot(dmAct)
%     drawnow
    toc
    
    tmp = uint8(tmp);
    clear tmp2
    
    outObj = VideoWriter([ms.dirName '/' 'msvideo.avi'],'Grayscale AVI');
    open(outObj)
    writeVideo(outObj,tmp);
    close(outObj)
    
end