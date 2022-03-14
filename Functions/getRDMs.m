function RDMs = getRDMs(envs,envLabel,morphPoint)
    features = [];
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        for k = 1:length(envI)
            features(envI(k),:) = [k i];
        end
    end
    lagMat = abs(bsxfun(@minus,[1:length(envs)],[1:length(envs)]'));
    envMat = abs(bsxfun(@minus,features(:,2),features(:,2)'));

    vec = [];
    for gi = 1:6
        vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
    end
    attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
        [features(:,2)'>=vec(1:length(features(:,2)))]'));

    RDMs = cat(3,lagMat,attractorMat,envMat);
    
end