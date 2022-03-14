function colors = factor3dcm(r2penalty)
    colors = bsxfun(@times,[1 0 0],r2penalty(:,1));
    colors = colors+bsxfun(@times,[0 0 1],r2penalty(:,2));
    colors = colors+bsxfun(@times,[0 1 0],r2penalty(:,3));
    colors(colors<0) = 0;
    colors = colors.^(1./4);
    colors = [colors./nanmax(colors(:))];
end