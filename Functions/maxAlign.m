function [maxA minA] = maxAlign(v,doMax)
    
    [blah bestSession] = nanmax(v,[],2);
    maxA = [nan(size(v)) v];
    for k = 1:length(v(:,1))
        maxA(k,:,:) = circshift(maxA(k,:),[0 -(bestSession(k)-1)]);
    end
    
    [blah bestSession] = nanmin(v,[],2);
    minA = [nan(size(v)) v];
    for k = 1:length(v(:,1))
        minA(k,:,:) = circshift(minA(k,:),[0 -(bestSession(k)-1)]);
    end
    
%     for k = 1:length(v(:,1))
%         [blah bestSession(k)] = find(~isnan(v(k,:)),1,'first');
%     end
%     firstA = [nan(size(v)) v];
%     for k = 1:length(v(:,1))
%         firstA(k,:,:) = circshift(firstA(k,:),[0 -(bestSession(k)-1)]);
%     end
%     
%     for k = 1:length(v(:,1))
%         [blah bestSession(k)] = find(~isnan(v(k,:)),1,'last');
%     end
%     lastA = [nan(size(v)) v];
%     for k = 1:length(v(:,1))
%         lastA(k,:,:) = circshift(lastA(k,:),[0 -(bestSession(k)-1)]);
%     end
end