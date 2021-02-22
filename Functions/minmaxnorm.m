function a = minmaxnorm(a)
    a = a - repmat(nanmin(nanmin(a,[],1),[],2),[size(a(:,:,1))]);
    norm = nanmax(nanmax(a,[],1),[],2);
    a = a./repmat(norm,[size(a(:,:,1))]);
end