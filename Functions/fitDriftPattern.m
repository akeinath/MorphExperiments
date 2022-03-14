function [params error] = fitDriftPattern(v,n)

    if nargin < 2 || isempty(n)
        n = 1;
    end

    nsims = 20;

    numParams = (n).*2+1;
    params = nan(length(v(1,1,:)),numParams);
    error = nan(length(v(1,1,:)),1);
    clearLen = 0;
    parfor i = 1:length(v(1,1,:))
%         str = sprintf('%0.3f',i.*100./length(v(1,1,:)));
%         fprintf([repmat('\b',[1 clearLen]) str]);
%         clearLen = length(str);
        
        x = [1:length(v(:,1,1))];
        
        pparams = nan(nsims,5);
        perror = nan(nsims,1);
        for si = 1:nsims
            seed = rand(1,numParams).*4-2;
            [pparams(si,:) perror(si)] = fmincon(@(p)DOL_error(v(:,:,i),p),seed ,[],[],[],[],...
                [ones(1,numParams).*-5],[ones(1,numParams).*5],[],optimoptions('fmincon','Display','none'));
        end
        [a b] = nanmin(perror);
        params(i,:) = pparams(b(1),:);
        error(i) = perror(b(1));
        
%         close all
%         figure(1)
%         set(gcf,'position',[50 50 400 200])
%         subplot(1,2,1)
%         imagesc(squarify(v(:,:,i)))
%         alpha(double(~isnan(squarify(v(:,:,i)))))
%         caxis([-1 1])
%         axis equal
%         axis off
%         subplot(1,2,2)
%         toPlot = DOL(x,params(i,:));
%         mask = triu(true(size(v(:,:,1))),1);
%         tv = toPlot(mask);
%         toPlot = toPlot';
%         toPlot(mask) = tv;
%         imagesc(toPlot)
%         alpha(double(~isnan(squarify(v(:,:,i)))))
%         caxis([-1 1])
%         axis equal
%         axis off
%         saveFig(gcf,['Plots/BatchedAnalyses/Summary/Clustering/DriftFitMatrices/' ...
%             num2str(i)],[{'pdf'} {'tiff'}]);
    end

end

function error = DOL_error(data,p)
    pred = DOL(1:length(data(1,:)),p);
    error = nanmean(((data(:)-pred(:)).^2));
end








