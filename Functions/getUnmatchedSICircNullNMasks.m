function [pvals actual null ] = getUnmatchedSICircNullNMasks(P,T,masks)
    pvals = nan;
    actual = nan;
    null = nan;
    velThresh = 2;
    nsims = 1000; %500
    minShift = 900;
    
    actual = help_getSI(P,T,masks);
    
    null = nan([length(actual) nsims]);
    parfor si = 1:nsims
        gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);
        null(:,si) = help_getSI(P,gT,masks);
    end
    
    pvals = 1-nanmean(bsxfun(@gt,actual,null),2);
    pvals(isnan(actual)) = nan;
end