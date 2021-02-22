function a = maxnorm(a)
    norm = nanmax(nanmax(a,[],1),[],2);
    a = a./repmat(norm,[size(a(:,:,1))]);
end