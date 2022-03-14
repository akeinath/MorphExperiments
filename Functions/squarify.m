function mat = squarify(mat)
    mat = nanmax(mat,permute(mat,[2 1 3]));
end