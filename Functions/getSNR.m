function getSNR(paths)
    clc
    close all
    drawnow
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    fprintf(['Computing signal-to-noise ratio:\n']);
    for p = paths'
        s = load(p{1});
        fprintf(['\t' num2str(p{1}) '\n']) 
        
        if isfield(s.processed,'snr')
            s.processed = rmfield(s.processed,'snr');
        end
        
        t = s.calcium.FiltTraces';
        t = s.calcium.FiltTraces(s.processed.validTraceFrames(:,1),:)';
        isGood = find(diff(t(1,:),[],2)~=0);
        for k = 1:length(t(:,1))
            t(k,isGood(1):isGood(end)) = linterp(isGood,t(k,isGood),isGood(1):isGood(end));
        end
        
        t = detrend(t')'; % norm to the median
        
        inlog = @(x)(log(x./(1-x)));
        

        nFrames = round(0.4./(1./30));
        snr = nan(length(t(:,1)),1);
        shsnr = nan(length(t(:,1)),2);
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing in serial: '])
        for k = 1:length(t(:,1))
            iter = iter+1;
            fprintf(repmat('\b',[1 strLength]));
            str = sprintf([ num2str(iter) ' of ' num2str(length(t(:,1)))]);
            fprintf(str);
            strLength = length(str);
                
            bsd = nanstd(t(k,t(k,:)<0))./sqrt(1-2./pi);
            zi = t(k,:)./bsd;
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
        fprintf('\n')
        s.processed.snr.whole = snr;
        s.processed.snr.splithalf = shsnr;
        save(p{1},'-struct','s','-v7.3');
    end
end