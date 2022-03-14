function binarizeTrace(paths)

    windowSize = 30.*30;

    clc
    fprintf(['\t\nNorming traces with 2nd-order autoregressive model.\n\n'])
    warning off all
    for p = paths'
        s = load(p{1});
        fprintf(['\t\t' num2str(p{1}) '\n\t\t\tModeling spike trains:      '])
        
        t = s.calcium.FiltTraces(s.processed.validTraceFrames(:,1),:)';
        isGood = find(diff(t(1,:),[],2)~=0);
        for k = 1:length(t(:,1))
            t(k,isGood(1):isGood(end)) = linterp(isGood,t(k,isGood),isGood(1):isGood(end));
        end
        
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing in serial: '])
        for k = 1:length(t(:,1))
            iter = iter+1;
            fprintf(repmat('\b',[1 strLength]));
            str = sprintf([ num2str(iter) ' of ' num2str(length(t(:,1)))]);
            fprintf(str);
            strLength = length(str);
            
            dt = [0 conv(diff(t(k,:)),fspecial('gauss',[1 60],5),'same')];

            bsd = nanstd(dt(t(k,:)<0))./sqrt(1-2./pi); %%% Z-score derivative based on estimated from below baseline portion
            zdt = dt./bsd;
            
            t(k,:) = zdt > 2.5;   
        end
        fprintf('\n')
        s.processed.trace = t;
        s.processed.traceModel = {'Smooth dZ(t) Rise Detection'};
        
        save(p{1},'-struct','s','-v7.3');
    end
end