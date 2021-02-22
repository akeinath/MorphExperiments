function [snr shsnr] = msSNR(ms)
        
    t = ms.FiltTraces';
    isGood = find(diff(t(1,:),[],2)~=0);
    for k = 1:length(t(:,1))
        t(k,isGood(1):isGood(end)) = linterp(isGood,t(k,isGood),isGood(1):isGood(end));
    end

    t = detrend(t')'; % norm to the median

    inlog = @(x)(log(x./(1-x)));

    snr = nan(length(t(:,1)),1);
    shsnr = nan(length(t(:,1)),2);
    parfor k = 1:length(t(:,1))
        bsd = nanstd(t(k,t(k,:)<0))./sqrt(1-2./pi);
        zi = t(k,:)./bsd;
        nFrames = round(0.4./(1./30));
        pif = normcdf(zi);

        shiftVals = nan(nFrames,length(pif));
        for shiftI = 0:nFrames-1
            shiftVals(shiftI+1,:) = circshift(pif,-shiftI,2);
        end
        pmin = nanmin(prod(shiftVals(:,1:end-nFrames),1)).^(1./nFrames);
        snr(k) = -inlog(pmin);

        half = 1:length(shiftVals(1,:)) < length(shiftVals(1,:))./2;
        shsnr(k,:) = -inlog([nanmin(prod(shiftVals(:,half),1)).^(1./nFrames) ...
            nanmin(prod(shiftVals(:,~half),1)).^(1./nFrames)]);
    end
end