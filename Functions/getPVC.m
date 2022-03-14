function pvc = getPVC(um)

    tic
    pvc = nan([length(um(1,1,1,:)) length(um(1,1,1,:))]);
    for si = 1:length(um(1,1,1,:))
        for sj = si+1:length(um(1,1,1,:))
            tmp1 = um(:,:,:,si);
            tmp2 = um(:,:,:,sj);
            
            badCells = permute(all(all(isnan(tmp1),1),2) | all(all(isnan(tmp2),1),2),[3 1 2]);
            
            if all(badCells)
                continue
            end
            
            tmp1(:,:,badCells) = [];
            tmp2(:,:,badCells) = [];
            
            tmp1 = reshape(tmp1,[numel(tmp1(:,:,1)) length(tmp1(1,1,:))]);
            tmp2 = reshape(tmp2,[numel(tmp2(:,:,1)) length(tmp2(1,1,:))]);
            
            isBadPixels = all(isnan(tmp1),2) | all(isnan(tmp2),2);
            tmp1(isBadPixels,:) = [];
            tmp2(isBadPixels,:) = [];
            
            pvc(si,sj) = corr(tmp1(:),tmp2(:));
        end
    end     
    tmp = toc;
end
