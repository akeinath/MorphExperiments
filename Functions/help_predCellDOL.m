function [out MAE] = help_predCellDOL(vals,testParams)
    if nargin < 2 || isempty(testParams)
        testParams = [0:3];
    end
    
    x = [1:length(vals(:,1,1))];
    out = nan([size(vals,1) size(vals,2) size(vals,3) length(testParams)]);
    for numParams = testParams
        for i = 1:length(vals(:,1,1))
            for j = i+1:length(vals(1,:,1))
                v = vals;
                v(i,j,:) = nan;

                [params error] = fitDriftPattern(v,numParams);
                for k = 1:length(vals(1,1,:))
                    pred = DOL(x,params(k,:));
                    out(i,j,k,numParams==testParams) = pred(i,j);
                end
            end
        end
    end
    
    MAE = nan(length(out(1,1,:,1)),length(out(1,1,1,:)));
    for i = 1:length(out(1,1,1,:))
        [a b] = mat2lag((vals-out(:,:,:,i)));
        MAE(:,i) = nanmean(b.^2,2);
    end
end