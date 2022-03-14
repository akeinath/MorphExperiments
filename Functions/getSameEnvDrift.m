function params = getSameEnvDrift(across)
    isSameEnv = cellfun(@strcmp,repmat(permute(across.RDMs.envs,[1 3 2]),[1 length(across.RDMs.envs(:,1)) 1]), ...
        repmat(permute(across.RDMs.envs,[3 1 2]),[length(across.RDMs.envs(:,1)) 1 1]));
    
    isSameSize = nansum(nansum(isSameEnv(:,:,1),1),2);
    
    sameEnvVals = reshape(across.RDMs.whole(isSameEnv),[isSameSize length(isSameEnv(1,1,:))]);
    
    lag = repmat(abs(bsxfun(@minus,[1:length(isSameEnv(:,1,1))],[1:length(isSameEnv(:,1,1))]')), ...
        [1 1 length(sameEnvVals(1,:))]);
    
    sameEnvLag = reshape(lag(isSameEnv),[isSameSize length(isSameEnv(1,1,:))]);
    
    polyN = 1;
    params = nan(length(sameEnvVals(1,:)),polyN+1);
    for k = 1:length(sameEnvVals(1,:))
        v = sameEnvVals(:,k);
        lv = sameEnvLag(:,k);
        
        isGood = ~isnan(v);
        
        if nansum(isGood) < 6
            continue
        end
        
        params(k,:) = polyfit(lv(isGood),v(isGood),polyN);
    end
end