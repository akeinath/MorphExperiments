function sfpXVal(SFPs,vals,scalar)
    if nargin < 3 || isempty(scalar)
        scalar = ones(length(vals(:,1)),1);
    end

    while length(vals(1,:))<3
        vals = [vals zeros(length(vals(:,1)),1)];
    end
    
    SFPs(:,:,any(isnan(vals),2)) = [];
    scalar(any(isnan(vals),2),:) = [];
    vals(any(isnan(vals),2),:) = [];
    
    colors = bsxfun(@times,[0 0.0 1],vals(:,2));
    colors = colors+bsxfun(@times,[1 0 0],vals(:,1));
    colors = colors+bsxfun(@times,[0 1 0],vals(:,3));
    
    nSFPs = repmat(permute(scalar,[3 2 1]),size(SFPs(:,:,1))).* ...
        [SFPs./repmat(nanmax(nanmax(SFPs,[],1),[],2),size(SFPs(:,:,1)))];
    t = bsxfun(@times,permute(nSFPs,[1 2 4 3]),permute(colors,[4 3 2 1]));
    image(nanmax(t,[],4))
end