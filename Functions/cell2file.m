function cell2file(vals,outP)
    fid = fopen(outP,'w');
    for i = 1:length(vals(:,1))
        for j = 1:length(vals(1,:))
            fprintf(fid,num2str(vals{i,j}));
            fprintf(fid,['\t']);
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end