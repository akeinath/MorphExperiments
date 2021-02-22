function si = help_getSI(P,T,masks)
    for k = 1:length(masks)
        [blah(:,:,:,k) s(:,:,k) m(:,:,:,k)] = mkTraceMaps(P,T,masks{k});
    end
    m = reshape(m,[numel(m(:,:,1)) length(m(1,1,:))]);
    s = s(:);
%     m(s<15,:) = nan;
%     s(s<15) = nan;
    normS = s./nansum(s(:));
%     normM = bsxfun(@rdivide,m,nanmean(m,1));
    normM = bsxfun(@rdivide,30.*m,30.*nanmean(T,2)');
    si = nansum(bsxfun(@times,normS,30.*m.*log2(normM)),1)';
end