function normTrace(paths)
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
        
% % %         t = detrend(t')'; % norm to the median
% % %         
% % %         iter = 0;
% % %         strLength = 0;
% % %         fprintf(['\t\tComputing in serial: '])
% % %         for k = 1:length(t(:,1))
% % %             iter = iter+1;
% % %             fprintf(repmat('\b',[1 strLength]));
% % %             str = sprintf([ num2str(iter) ' of ' num2str(length(t(:,1)))]);
% % %             fprintf(str);
% % %             strLength = length(str);
% % %             
% % %             dt = [0 conv(diff(t(k,:)),fspecial('gauss',[1 60],5),'same')];
% % %             
% % %             bsd = nanstd(dt(t(k,:)<0))./sqrt(1-2./pi); %%% Z-score derivative based on estimated from below baseline portion
% % %             zdt = dt./bsd;
% % %             t(k,:) = zdt > 2.5;            
% % %         end
% % %         fprintf('\n')
% % %         s.processed.trace = t;
% % %         s.processed.traceModel = {'Smooth dZ(t) Rise Detection'};

        
%         t = detrend(t')'; % norm to the median
%         
%         inlog = @(x)(log(x./(1-x)));
%         
%         threshold = 0.25;
%         iter = 0;
%         strLength = 0;
%         fprintf(['\t\tComputing in serial: '])
%         for k = 1:length(t(:,1))
%             iter = iter+1;
%             fprintf(repmat('\b',[1 strLength]));
%             str = sprintf([ num2str(iter) ' of ' num2str(length(t(:,1)))]);
%             fprintf(str);
%             strLength = length(str);
%                 
%             bsd = nanstd(t(k,t(k,:)<0))./sqrt(1-2./pi);
%             zi = t(k,:)./bsd;
% %             plot(zi)
% %             hold on
% %             plot(diff(zi))
% %             plot(conv(diff(zi),fspecial('gauss',[1 60],5),'same'))
%             t(k,:) = [conv(diff(zi),fspecial('gauss',[1 60],5),'same')>threshold false];
%         end
%         fprintf('\n')
%         s.processed.trace = t;
%         s.processed.traceModel = {'Smooth dZ(t) Rise Detection'};

        
        
        t = t';
        spikes = nan(size(s.calcium.FiltTraces(s.processed.validTraceFrames(:,1),:)));
        for i = 1:length(s.calcium.FiltTraces(1,:))
            str_f = sprintf('%6.1f',100.*i/length(s.calcium.FiltTraces(1,:)));
            fprintf([repmat('\b',[1 6]) str_f])
            [blah spikes(:,i)] = deconvolveCa(t(:,i),'ar2');
        end
        
        fprintf('\n')
        s.processed.trace = spikes';
        s.processed.traceModel = {'AutoReg2'};
        
        save(p{1},'-struct','s','-v7.3');
    end
end