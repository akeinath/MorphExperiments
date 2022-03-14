function c = hulkcm(vals,colorLims)
    if nargin < 2 || isempty(colorLims)
        colorLims = [nanmin(vals(:)) nanmax(vals(:))];
    end
    
    nvals = (vals-colorLims(1))./colorLims(2);

    c = bsxfun(@times,nvals(:,1),[0.65 0.0 0.65])+0.1 + ...
            bsxfun(@times,nvals(:,2),[0 0.65 0])+0.1;
end