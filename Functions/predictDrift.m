function predictDrift(vals,root)

    testParams = [0:3];
    [out MAE] = help_predCellDOL(vals,testParams);
    [out MAE2] = help_predCellDOS(vals,[1 2]);
    
% % % %     mk Shuffled values within cell % shuffle rows/cols together randomly
    
    sVals = vals;
    for k = 1:length(vals(1,1,:))
        o = randperm(length(vals(:,1,1)));
        tmp = squarify(vals(o,o,k));
        tmp(isnan(vals(:,:,1))) = nan;
        sVals(:,:,k) = tmp;
    end
    [sout sMAE] = help_predCellDOL(sVals,testParams);
    [sout sMAE2] = help_predCellDOS(sVals,[1 2]);
%     
%     save('ContextPredictionActual','ae');
    
    ae = [];
    MAE = [MAE MAE2];
    sMAE = [sMAE sMAE2];
    for i = 1:length(MAE(1,:))
        ae = [ae MAE(:,i) sMAE(:,i)];
    end
    
    pairwiseNonparametrics(ae,[root '/MAE_']);
    
    close all
    figure(1)
    set(gcf,'position',[50 50 300 225])
    mkBow(log10(ae))
    ylabel('log_1_0(MAE)')
    xlabel('Order of DOL fit')
    saveFig(gcf,[root '/MAE_log'],[{'tiff'} {'pdf'}])    
    
    doPerms = perms(1:length(vals(:,1,1)));
    if exist('DriftPredictionNull.mat','file')==2
        load('DriftPredictionNull');
    else
        null_MAE = nan([length(vals(1,1,:)) length(testParams) ...
            length(doPerms(:,1))]);
    end
    fprintf(['\n\t\t\tFitting drift pattern for all perms:  '])
    clearLen = 0;
    for o = fliplr(doPerms')
        clc
        ind = find(ismember(doPerms,o','rows'));
        if ~isnan(null_MAE(1,1,ind))
            continue
        end
        str = sprintf([num2str(ind) ' of ' num2str(length(doPerms(:,1)))]);
        fprintf([repmat('\b',[1 clearLen]) str])
        clearLen = length(str);
        
        sVals = squarify(vals(o,o,:));
        sVals(isnan(vals)) = nan;
        tic
        [out sMAE] = help_predCellDOL(sVals,testParams);
        toc
        null_MAE(:,:,ind) = sMAE;
        save('DriftPredictionNull','null_MAE');
    end
    
    pVal = nanmean(bsxfun(@gt,null_MAE(:,:,end),null_MAE(:,:,1:end-1)),3);
    
    tmp = permute(nanmin(null_MAE(:,:,:),[],2),[1 3 2]);    
    cumHist(nanmean(bsxfun(@gt,tmp(:,end),tmp(:,1:end-1)),2),[0:0.01:1])
    
%     sVals = vals;
%     for k = 1:length(vals(1,1,:))
%         o = randperm(length(vals(:,1,1)));
%         tmp = squarify(vals(o,o,k));
%         tmp(isnan(vals(:,:,1))) = nan;
%         sVals(:,:,k) = tmp;
%     end
    [sout sMAE] = help_predCellDOL(sVals,testParams);
    
    
    
    ae = [];
    for i = 1:length(testParams)
        ae = [ae MAE(:,i) sMAE(:,i)];
    end
   
    pairwiseNonparametrics(ae,[root '/MAE_'])
    
    close all
    figure(1)
    set(gcf,'position',[50 50 300 225])
    mkBow(log10(ae))
    ylabel('log_1_0(MAE)')
    xlabel('Order of DOL fit')
    saveFig(gcf,[root '/MAE_log'],[{'tiff'} {'pdf'}])
    
%     %%%%% Surrogate data

    for surrogateParams = testParams
        [params error] = fitDriftPattern(vals,surrogateParams);
        fakeV = nan(size(vals));
        x = [1:length(vals(:,1,1))];
        for k = 1:length(vals(1,1,:))
            fakeV(:,:,k) = DOL(x,params(k,:));
        end
        fakeV(isnan(vals)) = nan;
        fakeV = fakeV+randn(size(fakeV)).*0.05; % add a little noise

        sfVals = fakeV;
        for k = 1:length(vals(1,1,:))
            o = randperm(length(fakeV(:,1,1)));
            tmp = squarify(fakeV(o,o,k));
            tmp(isnan(fakeV(:,:,1))) = nan;
            sfVals(:,:,k) = tmp;
        end

        [out MAE] = help_predCellDOL(fakeV,testParams);
        [sout sMAE] = help_predCellDOL(sfVals,testParams);
        
        ae = [];
        for i = 1:length(testParams)
            ae = [ae MAE(:,i) sMAE(:,i)];
        end
    
        close all
        figure(1)
        set(gcf,'position',[50 50 300 225])
        mkBow(log10(ae))
        set(gca,'ylim',[-2 0.5])
        ylabel('log_1_0(MAE)')
        xlabel('Order of DOL fit')
        saveFig(gcf,[root '/SurrogateData/MAE_log_SurrogateParams_' ...
            num2str(surrogateParams)],[{'tiff'} {'pdf'}])
    end
end











































