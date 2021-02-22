function classXDynamics(d,doFold)

    hold on
    colors = [0.9 0.2 0.2; 0.2 0.2 0.9; 0.2 0.7 0.2];
    for k = 1:length(d)
        
        if length(d{k}(1,:)) < 64
            d{k} = [nan(length(d{k}(:,1)),32) d{k}];
        end
        
        if nargin > 1 && doFold
            set(gca,'xlim',[0 31])
            d{k} = [nan(length(d{k}(:,1)),32) nanmean(cat(3,d{k}(:,33:end),fliplr(d{k}(:,1:32))),3)];
        end
        
%         qrts = nan(5,length(d{k}(1,:)));
%         for j = 1:length(d{k}(1,:))
%             tmp = sort(d{k}(~isnan(d{k}(:,j)),j));
%             if ~isempty(tmp)
%                 qrts(:,j) = tmp([1 round(length(tmp).*[0.25:0.25:0.75]) end]);
%             end
%         end
        
        x = [-32:31];
        
        m = nanmean(d{k});
        se = nanstd(d{k})./sqrt(nansum(~isnan(d{k})));
        patch([x(1:end-1); x(1:end-1); x(2:end); x(2:end)], ...
            [m(1:end-1)-se(1:end-1); m(1:end-1)+se(1:end-1); ...
            m(2:end)+se(2:end); m(2:end)-se(2:end)],colors(k,:),...
            'edgecolor','none','facealpha',0.25)

        plot(x,m+se,'color',colors(k,:).*0.5+0.5);
        plot(x,m-se,'color',colors(k,:).*0.5+0.5);
        plot(x,m,'color',colors(k,:));

        
%         hold on
%         plot(x,qrts(2,:),'color',colors(k,:).*0.5+0.5);
%         plot(x,qrts(4,:),'color',colors(k,:).*0.5+0.5);
%         plot(x,qrts(3,:),'color',colors(k,:));
%         plot(x,qrts(1,:),'color',colors(k,:).*0.25+0.75);
%         plot(x,qrts(5,:),'color',colors(k,:).*0.25+0.75);
        
        
%         patch([x(1:end-1); x(1:end-1); x(2:end); x(2:end)], ...
%             [m(1:end-1)-se(1:end-1); m(1:end-1)+se(1:end-1); ...
%             m(2:end)+se(2:end); m(2:end)-se(2:end)],colors(k,:),...
%             'edgecolor','none','facealpha',0.5)
    end
    

% % %     if nargin > 1 && doFold
% % % %         for k = 1:length(d)
% % % %             tmpV = d{k}(:,34:end);
% % % %             tmpX = repmat(x(34:end),[length(d{k(:,1)}) 1]);
% % % %             [rval pval] = corr(tmpX(~isnan(tmpV)),...
% % % %                 tmpV(~isnan(tmpV)));
% % % % 
% % % %             textX = get(gca,'xlim');
% % % %             textX = textX(2) - [textX(2)-textX(1)].*0.9;
% % % %             textY = get(gca,'ylim');
% % % %             textY = textY(1) + [textY(2)-textY(1)].*0.1.*(4-k);
% % % %             text(textX,textY,sprintf(['r = %.3f, p = %.2e'],[rval pval]), ...
% % % %                 'fontweight','normal','fontname','arial','fontsize',10);
% % % %         end
% % % 
% % %          tmpX = x(1,34:end);
% % %          comb = 0;
% % %          for j = 1:length(d)
% % %               for k = j+1:length(d)
% % %                   a = d{j}(:,34:end);
% % %                   b = d{k}(:,34:end);
% % %                   pvals = nan(1,length(a(1,:)));
% % %                   for i = 1:length(pvals(1,:))
% % %                       if ~isempty(a(~isnan(a(:,i)),i)) && ~isempty(b(~isnan(b(:,i)),i))
% % %                         pvals(i) = ranksum(a(~isnan(a(:,i)),i),b(~isnan(b(:,i)),i));
% % %                       end
% % %                   end
% % % 
% % %                   if any(pvals<0.01)
% % %                       textY = get(gca,'ylim');
% % %                       textY = textY(2) - [textY(2)-textY(1)].*0.05 - [textY(2)-textY(1)].*0.05.*(comb);
% % % 
% % %                       plot(tmpX(pvals<0.05),textY,'linestyle','none','marker','o', ...
% % %                           'color',colors(j,:),'markerfacecolor',colors(k,:),'linewidth',2);
% % %                       comb = comb+1;
% % %                   end
% % %               end
% % %          end
% % %     end
end


















