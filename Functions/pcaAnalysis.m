function pcaAnalysis(sim,envs,envLabel)
    
    goodPix = false(size(sim(:,:,1)));
    goodPix = triu(~goodPix,1);
    
    rs = reshape(sim(repmat(goodPix,[1 1 length(sim(1,1,:))])), ...
        [nansum(goodPix(:)) length(sim(1,1,:))]);
    
    [a b c d e] = pca(rs');
    
    pcMaps = nan(size(sim));
    pcMaps(repmat(goodPix,[1 1 length(sim(1,1,:))])) = b';
    
    mds2D(squarify(pcMaps(:,:,4)),envs,envLabel)
end