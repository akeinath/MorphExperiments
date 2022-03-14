function RDMs = getRDMs_r2(envs,envLabel,varargin)
    features = [];
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        for k = 1:length(envI)
            features(envI(k),:) = [k i];
        end
    end
    lagMat = abs(bsxfun(@minus,[1:length(envs)],[1:length(envs)]'));
    envMat = abs(bsxfun(@minus,features(:,2),features(:,2)'));

    allAttractorMats = [];
    for morphPoint = 2:6
        vec = [];
        for gi = 1:6
            vec((gi-1).*6+1:(gi).*6) = morphPoint;
        end
        attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
            [features(:,2)'>=vec(1:length(features(:,2)))]'));
        allAttractorMats = cat(3,allAttractorMats,attractorMat);
    end

    RDMs = cat(3,lagMat,lagMat.^2,lagMat.^4,allAttractorMats);
end