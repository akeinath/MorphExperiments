function [ap] = fitSigmoidTimeDecay(sim)

%     sim = atanh(sim);

    rsim = [nan(size(sim)) sim];
    for i = 1:length(sim(:,:,1))
        rsim(i,:,:) = circshift(rsim(i,:,:),[0 -(i-1)]);
    end
    val = permute(nanmean(rsim(:,33:end,:),1),[3 2 1]);
    
     
    ap = nan(length(val(:,1)),3);
    x = [1:length(val(1,:))]-1;
    for k = 1:length(val(:,1))
        y = val(k,:);
        
        [params fval] = fmincon(@(p)sig_error(y,p),[(x(end)./2)+0.5 1 0 2] ,[],[],[],[],...
            [-inf 1e-10 -1 0],[inf 100 1 1],[],optimoptions('fmincon','Display','none'));
        
        tmpA = doSig(0:31,params);
        
        lA = polyval([diff(tmpA([1 end]))./length(tmpA) tmpA(1)],[0:length(tmpA)-1]);

        nTmpA = (tmpA-tmpA(end))./range(tmpA);
        nLA = (lA-lA(end))./range(lA);

        ap(k,1) = nanmean(abs([nTmpA-nLA])); % mean squared difference
        
        ap(k,2) = 1 - nansum((y-tmpA).^2)./nansum((x-nanmean(x)).^2);
        ap(k,3) = 1 - nansum((y-lA).^2)./nansum((x-nanmean(x)).^2);
    end
    
    
%     ap = nan(length(rsim(1,1,:)),2);
%     x = [1:length(rsim(1,:,1))]-1;
%     for k = 1 %:length(rsim(1,1,:))
%         y = rsim(:,:,k);
%         
%         
%         [params fval] = fmincon(@(p)sig_time_error(y,p),[(x(end)./2)+0.5 1 0 2] ,[],[],[],[],...
%             [-inf 1e-10 -1 0],[inf 10 1 1],[],optimoptions('fmincon','Display','none'));
%         
%         tmp = doSigTime(0:31,params)
%         
%         [p S] = polyfit(x(~isnan(y)),y(~isnan(y)),1);
%         ap(k,:) = fliplr(p);
%     end
end

function error = sig_error(data,p)
    pred = doSig(1:length(data(1,:)),p);
    error = nanmean((data-pred).^2);
end

function pred = doSig(x,p)
    pred = ([1./(exp((-fliplr(x)+p(1)).*p(2))+1)].*p(4))+p(3);
end