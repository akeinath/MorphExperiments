function conXdrift(sim,um,envs,root,mfr)
    

    fprintf(['\n\n\t Aggregated context vs drift analysis...  '])
    strLen = 0;
    
    conMaps = [];
    catDat = [];
    aCon = [];
    for doConI = 1:5
        str = sprintf('Context Sequence %i',doConI);
        fprintf([repmat('\b',[1 strLen]) str])
        strLen = length(str);
        
        %%%% Compute lagged drift, excluding the context sequence
        rawDrift = nan(5,5,size(sim,5));
        doDrift = setxor([1:5],doConI);
        for doDriftI = 1:length(doDrift)
            for doDriftJ = doDriftI+1:length(doDrift)
                
                tmp = permute(sim(:,:,doDrift(doDriftI),doDrift(doDriftJ),:),[1 2 5 3 4]);
                
                rt = reshape(tmp(repmat(logical(eye(6)),[1 1 length(tmp(1,1,:))])), ...
                    [6 length(tmp(1,1,:))]);
                rawDrift(doDrift(doDriftI),doDrift(doDriftJ),:) = nanmean(rt,1);
            end            
        end
        
        lagDrift = nan(size(rawDrift,3),4);
        for k = 1:size(rawDrift,3)
            tmp = rawDrift(:,:,k);
            for i = 1:4
                lagDrift(k,i) = nanmean(diag(tmp,i));
            end
        end
        
        %%% Compute context sequence fit
        
        conSim = permute(sim(:,:,doConI,doConI,:),[1 2 5 3 4]);
        aCon = cat(3,aCon,conSim);
        
        conMaps = [cat(3,conMaps,um(:,:,:,(doConI-1).*6+1:(doConI).*6))];
        goodPixels = logical(triu(ones(6),1));
        
        rc = reshape(conSim(repmat(goodPixels,[1 1 size(conSim,3)])), ...
            [nansum(goodPixels(:)) size(conSim,3)]);
        doCon = ~any(isnan(rc),1);
        
        [contextParams contextError] = fitContextPattern(conSim(:,:,doCon),1); % og 50 sims
        conMSE = nan(size(conSim,3),1);
        conMSE(doCon) = contextError;   
        conPar = nan(size(conSim,3),5);
        conPar(doCon,:) = contextParams;
        catDat = [catDat; conMSE lagDrift conPar];
    end
    
%     save('catDat','catDat');
    
%     load('catDat')

    figure
    set(gcf,'position',[50 50 1200 250])
    for i = 1:4
        isGood = ~isnan(catDat(:,1)) & ~isnan(catDat(:,i+1));
        
        subplot(1,4,i)
        scatter(-log10(catDat(isGood,1)),catDat(isGood,i+1),3,[0.5 0.5 0.5]);
        lsline
        set(gca,'ylim',[-1 1],'xlim',[0 5])
        hold on
        plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--');
        [r p] = corr(-catDat(isGood,1),catDat(isGood,i+1),'type','spearman');
        text(4.75,-0.7,sprintf('rho = %0.3f\np = %0.3e\nn = %i',[r p nansum(isGood)]),...
            'fontname','arial','fontsize',9,'horizontalalignment','right')
        axis square
        xlabel('-MSEContext')
        ylabel('Stability (r)')
    end
    saveFig(gcf,root,[{'pdf'} {'tiff'}]);
    
    %%%%%% Partition MFR-matching analysis
    
    if nargin > 2 && ~isempty(mfr)
        
        nParts = 6;
        
        groupMFR = [];
        for i = 1:5
            doS = (i-1).*6+1:(i).*6;
            groupMFR = [groupMFR; nanmean(mfr(:,doS),2)];
        end
        
        matchedLagDat = [];
        
        av = [];
        for i = 1:5
            [rval pval] = corr(groupMFR(~isnan(catDat(:,1))&~isnan(catDat(:,i)),1), ...
                catDat(~isnan(catDat(:,1))&~isnan(catDat(:,i)),i),'type','spearman');
            av = [av; [rval pval nansum(~isnan(catDat(:,1))&~isnan(catDat(:,i)))]];
        end
        fprintf(sprintf('\nFiring rate covariates: \n'));
        fprintf(sprintf('\nr = %0.3f, p = %0.3e, n = %i',av'));
        fprintf(sprintf('\n'));
        
        for j = 1:4
            isGood = ~isnan(catDat(:,1)) & ~isnan(catDat(:,j+1));
        
            goodCD = catDat(isGood,[1 j+1]);
            goodMFR = groupMFR(isGood);
            
            [a b] = sort(goodCD(:,1),'ascend');
            cutoffs = [ones(1,nParts) 0]+length(b).*([0:nParts]./(nParts));
            groupedInds = [];
            groupedMFR = [];
            for i = 1:nParts
                groupedInds = [groupedInds {b(cutoffs(i):cutoffs(i+1))}];
                groupedMFR = [groupedMFR {goodMFR(b(cutoffs(i):cutoffs(i+1)))}];
            end
            
            figure
            set(gcf,'position',[50 50 600 200])
            subplot(1,2,1)
            [blah matchedInds] = compHist(groupedMFR,[0:0.0025:0.08],1);
            set(gca,'ylim',[0 100])
            
            nmi = nan(size(matchedInds));
            aggDat = [];
            for i = 1:nParts
                nmi(:,i) = groupedInds{i}(matchedInds(:,i));
                aggDat = [aggDat; {goodCD(nmi(:,i),2)}];
            end
            
            matchedLagDat = [matchedLagDat aggDat];
            
            tmp = [];
            for i = 1:nParts
                tmp = [tmp {goodMFR(nmi(:,i))}];
            end
            subplot(1,2,2)
            compHist(tmp,[0:0.0025:0.08],1);
            set(gca,'ylim',[0 100])
            saveFig(gcf,[root '_matchingExample'],[{'tiff'} {'pdf'}])
        end
        
        figure
        set(gcf,'position',[50 50 550 250])
        mkBar(matchedLagDat)
        set(gca,'ylim',[-1 1])
        hold on
        plot(get(gca,'xlim'),[0 0],'color','k','linestyle','--')
        xlabel('Lag')
        ylabel('Stability (r)')
        saveFig(gcf,[root '_matchedMFR'],[{'tiff'} {'pdf'}])
        
        cellOut = [];
        fid = fopen([root '_stats.txt'],'w');
        fprintf(fid,'\n\n\t\t%%%%%%%%%% Matched MFR Grouped-Con, Stability Analysis %%%%%%%%%%\n');
        for q = 1:4
            cellOut = [cellOut; repmat({[]},[1 4])];
            fprintf(fid,['\n\t\t\t Lag ' num2str(q)]);
            [pval,tbl] = anova1(cat(2,matchedLagDat{:,q}),[],'off');
            fprintf(fid,'\n\t F(%i,%i) = %0.2f, p = %.2e',[tbl{2,3}; tbl{3,3}; tbl{2,5}; tbl{2,6}]);
            tmp = cat(2,matchedLagDat{:,q});
            cellOut = [cellOut; [{sprintf('%i',q)} {'ANOVA'} {'all'} ...
                {sprintf('F(%i,%i) = %0.2f, p = %.2e',[tbl{2,3}; tbl{3,3}; tbl{2,5}; tbl{2,6}])}]];
            for i = 1:length(tmp(1,:))
                for j = i+1:length(tmp(1,:))
                    [a b c d] = ttest2(tmp(:,i),tmp(:,j));
                    fprintf(fid,'\n\t Comparison(%i, %i) t(%i) = %0.2f, p = %.2e',[i; j; d.df; d.tstat; b]);
                    cellOut = [cellOut; [{sprintf('%i',q)} {'t-test'} ...
                        {sprintf('%i vs. %i',[i; j])} ...
                        {sprintf('t(%i) = %0.2f, p = %.2e',[d.df; d.tstat; b])}]];
                end
            end
        end
        fclose(fid);
        
        if exist([root '_stats.xls'])==2
            delete([root '_stats.xls'])
        end
        xlswrite([root '_stats.xls'],cellOut);
    end
    
    %%% RSM dim reduction
    
% % %     load('ParExample_RSM')
% % %     
% % %     parDat = catDat(:,end-4:end);
% % %     gpd = parDat(~any(isnan(parDat),2),:);
% % %     gpCon = aCon(:,:,~any(isnan(parDat),2));
% % %     gpE = catDat(~any(isnan(parDat),2),1);
% % %     gpUm = conMaps(:,:,~any(isnan(parDat),2),:);
% % %     
% % %     [a b] = sort(gpE);
% % %     isGood = b(1:round(length(b).*1)); %%%% Best XX% of fits
% % %     
% % %     gpd = gpd(isGood,:);
% % %     gpCon = gpCon(:,:,isGood);
% % %     gpE = gpE(isGood,:);
% % %     gpUm = gpUm(:,:,isGood,:);
% % %     
% % %     out = clustParams(gpd,1:length(gpd),'mahalanobis',12); % 'minkowski'
% % %     
% % %     doColors = magma(length(out.inds)+1);
% % %     
% % %     figure
% % %     set(gcf,'position',[50 50 500 500])
% % %     for k = 1:length(out.inds)
% % %         plot(out.embedding(out.inds{k},1),out.embedding(out.inds{k},2), ...
% % %             'linestyle','none','marker','o','color',doColors(k,:),'markersize',3);
% % %         hold on
% % %     end
% % %     
% % %     figure
% % %     set(gcf,'position',[50 50 ceil(sqrt(length(out.inds))).*150 ceil(sqrt(length(out.inds))).*150])
% % %     for k = 1:length(out.inds)
% % %         subplot(ceil(sqrt(length(out.inds))),ceil(sqrt(length(out.inds))),k)
% % %         tmp = squarify(nanmean(gpCon(:,:,out.inds{k}),3));
% % %         imagesc(tmp)
% % %         caxis([0 1])
% % %         colormap inferno
% % %         alpha(double(~isnan(tmp)))
% % %     end
% % %     
% % %     doCLims = [[0 0 -0.5 0 -5]',[7 10 0.5 10 5]'];
% % %     figure
% % %     set(gcf,'position',[50 50 1500 600])
% % %     tl = [{'Transition Point (p1)'} {'Transition Abruptness (p2)'} ...
% % %         {'Assymetry (p3)'} {'Degree of modulation (p4)'} ...
% % %         {'Average overall similarity (p5)'}];
% % %     ctl = [0 7; 0 10; -0.5 0.5; 0 10; -5 5];
% % %     lim = round(range(out.embedding)./2).*1.1;
% % %     for i = 1:5
% % %         subplot(2,3,i)
% % %         doColors = v2rgb(gpd(:,i),ctl(i,:));
% % %         scatter(out.embedding(:,1), ...
% % %             out.embedding(:,2),10,doColors,'filled');
% % % %         set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
% % %         title(tl{i})
% % %         caxis(ctl(i,:))
% % %         axis square
% % %         axis equal
% % %     end
% % %     subplot(2,3,6)
% % %     plotEmbeddingExamples(out,ex)
% % %     title('Examples')
% % % %     set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
% % %     axis square
% % %     axis equal
% % % %     doColors = v2rgb(-log10(gpE));
% % % %     scatter(out.embedding(:,1), ...
% % % %         out.embedding(:,2),10,doColors,'filled');
% % % %     lim = (nanmax(abs([get(gca,'xlim') get(gca,'ylim')]))).*0.8;
% % % %     set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
% % % %     axis square
% % %     saveFig(gcf,[root '/Drift_UMAP_ParameterProjection'],[{'pdf'} {'tiff'}])
% % %     
% % % %     getEmbeddingExamples(out,gpCon);
% % % 
% % % %     doK = [8 4];
% % % %     for pi = 1:floor(length(gpCon(1,1,:))./prod(doK))
% % % %         tmp = gpCon(:,:,(pi-1).*prod(doK)+1:nanmin((pi).*prod(doK),length(gpCon(1,1,:))));
% % % %         
% % % %         close all
% % % %         drawnow
% % % %         figure
% % % %         set(gcf,'position',[50 50 800 600])
% % % %         for k = 1:length(tmp(1,1,:))
% % % %             subplot(4,8,k)
% % % %             imagesc(squarify(tmp(:,:,k)))
% % % %             alpha(double(~isnan(squarify(tmp(:,:,k)))))
% % % %             caxis([-1 1])
% % % %             axis off
% % % %             title(num2str((pi-1).*prod(doK)+k),'fontname','arial','fontsize',9)
% % % %             colormap magma
% % % %             axis square
% % % %             axis equal
% % % %         end
% % % %         outP = [root '/ContextualRSMs_' num2str(pi)];
% % % %         saveFig(gcf,outP,[{'pdf'} {'tiff'}]);
% % % %     end
% % %     
% % % %     doK = [2 4];
% % % %     for pi = 1:floor(length(gpCon(1,1,:))./prod(doK))
% % % %         tmp = gpUm(:,:,(pi-1).*prod(doK)+1:nanmin((pi).*prod(doK),length(gpCon(1,1,:))),:);
% % % %         
% % % %         close all
% % % %         drawnow
% % % %         figure
% % % %         set(gcf,'position',[50 50 800 600])
% % % %         for k = 1:length(tmp(1,1,:,1))
% % % %             catTmp = permute(tmp(:,:,k,:),[1 2 4 3]);
% % % %             catTmp = maxnorm(catTmp);
% % % %             blah = [];
% % % %             for q = 1:length(catTmp(1,1,:))
% % % %                 blah = [blah catTmp(:,:,q) nan(length(catTmp(:,1,1)),1)];
% % % %             end
% % % %             catTmp = blah(:,1:end-1);
% % % %             
% % % %             subplot(4,2,k)
% % % %             imagesc(catTmp)
% % % %             alpha(double(~isnan(catTmp)))
% % % %             caxis([0 1])
% % % %             axis off
% % % %             title(num2str((pi-1).*prod(doK)+k),'fontname','arial','fontsize',9)
% % % %             colormap parula
% % % %             axis square
% % % %             axis equal
% % % %         end
% % % %         outP = [root '/ContextualRSMs_' num2str(pi)];
% % % %         saveFig(gcf,outP,[{'pdf'} {'tiff'}]);
% % % %     end
% % %     
% % %     
% % % %     chooseEmbeddingExamples(out,gpCon);
end
































