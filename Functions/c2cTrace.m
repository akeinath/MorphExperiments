function c2cTrace(uGT)

    for si = 1:length(uGT)
        t = uGT{si};
        smth = 30;
        kern = fspecial('gauss',[1 smth.*7],smth);
        st = imfilter(t,kern,'same');
        
        goodCells = ~any(isnan(st),2);
        
        st(~goodCells,:) = [];
        
        pop = nanmean(st);
        
        c2c = nan(length(st(:,1)),length(st(:,1)));
        for ki = 1:length(st(:,1))       
            tic
            for kj = ki+1:length(st(:,1))
                pop = nanmean(st([1:ki-1 ki+1:end],:));
                beta = glmfit([st(ki,:); pop]',st(kj,:)');
                c2c(ki,kj) = beta(2);
                
                ap = c2c(~isnan(c2c));                
            end
            toc
            figure
            set(gcf,'position',[50 50 300 225])
            histogram(ap,[-0.2:0.025:0.5])
            set(gca,'xlim',[-0.2 0.5])
            drawnow
        end
    end

end