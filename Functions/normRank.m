function nv = normRank(v,dir)

    if nargin < 2 || isempty(dir)
        dir = 'descend';
    end
    
    nv = nan(size(v));
    [a sortInds] = sort(v,dir);
    for i = 1:length(v(1,:))
        [a nv(:,i)] = ismember(1:length(v(:,1)),sortInds(:,i));
    end
end