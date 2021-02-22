function warpMap = pvWarp(a,b)
    warpMap = [];
    
    % eliminate unregistered cells
    goodCells = ~all(all(isnan(a),1),2) & ~all(all(isnan(b),1),2);
    a(:,:,~goodCells) = [];
    b(:,:,~goodCells) = [];
    
    ra = reshape(a,[numel(a(:,:,1)) length(a(1,1,:))]);
    rb = reshape(b,[numel(b(:,:,1)) length(b(1,1,:))]);
    
    pix2pix = corr(ra',rb');
    
    [val indA] = nanmax(pix2pix,[],2);
    
    linMapA = [[1:length(ra(:,1))]' indA];
    
    goodInds = ~all(isnan(ra),2);
    [xa1 ya1] = ind2sub(size(a(:,:,1)),linMapA(goodInds,1));
    [xa2 ya2] = ind2sub(size(a(:,:,1)),linMapA(goodInds,2));
    
    doC = v2rgb(val(goodInds),[-1 1]);
    
    for k = 1:length(xa1)
        hold on
        plot([xa1(k) xa2(k)]',[ya1(k) ya2(k)]','color',doC(k,:),'marker','none',...
            'markerfacecolor',doC(k,:),'markersize',5);
    end
    set(gca,'xlim',[0 17],'ylim',[0 17])
end