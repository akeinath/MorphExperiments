function out = normCrop(v)
    normer = repmat(permute(nanmax(v,[],2),[2 3 1]),[length(v(1,:)) length(v(1,:))]);
    out = 1-[abs(bsxfun(@minus,permute(v,[3 2 1]),permute(v,[2 3 1])))./normer];
    out(repmat(tril(true(size(out(:,:,1)))),[1 1 length(out(1,1,:))])) = nan;
end