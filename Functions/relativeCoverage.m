function rc = relativeCoverage(aSamp,minSamp,envs)
    rc = nan(1,length(aSamp(1,1,:)));
    for i = 1:32
        rc(i) = nansum(nansum(aSamp(:,:,i)>=minSamp,1),2) ./ ...
            nansum(nansum(any(aSamp(:,:,ismember(envs,envs(i)))>0,3),1),2); % Half second sampling
    end
end