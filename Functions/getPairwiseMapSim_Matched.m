function out = getPairwiseMapSim_Matched(aSamp,um,uP,uGT,varargin)
    out = struct;
    
    if ismember({'pearson'},varargin)
        %%% Pearson correlation similarity
        os = round(aSamp.*30);
        nsims = 1;
        strLength = 0;
        
        tmp = (~permute(all(all(isnan(um),1),2),[3 4 1 2]));
        counts = nan(length(tmp(1,:)));
        for si = 1:length(tmp(1,:))
            for sj = si+1:length(tmp(1,:))
                counts(si,sj) = nansum(tmp(:,si)&tmp(:,sj));
            end
        end
        isGood = triu(true(size(counts)),1);
        matchNum = nanmin(counts(isGood));
        
        fprintf('\n\tComputing map similarity (matched sampling)...  ')
        tic
        pearson = nan([length(uGT) length(uGT) length(uGT{1}(:,1,1))]);
        pv = nan([length(uGT) length(uGT)]);
        zeroedpv = nan([length(uGT) length(uGT)]);
        for si = 1:length(uGT)
            for sj = si+1:length(uGT)
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf(['(' num2str(si) ', ' num2str(sj) ') of ('  ...
                    num2str(length(uGT)) ', ' num2str(length(uGT)) ')']);
                fprintf(str);
                strLength = length(str);
                
                tos = nanmin(os(:,:,si),os(:,:,sj));
                
                m1 = mkTraceMaps(uP{si},uGT{si},[],size(tos),tos,nsims);
                m2 = mkTraceMaps(uP{sj},uGT{sj},[],size(tos),tos,nsims);
                
                rt1 = reshape(m1,[numel(m1(:,:,1,1)) length(m1(1,1,:,1))]);
                rt2 = reshape(m2,[numel(m2(:,:,1,1)) length(m2(1,1,:,1))]);
                isGood = ~all(isnan(rt1),2) & ~all(isnan(rt2),2); 
                
                art1 = rt1(isGood,:);
                art2 = rt2(isGood,:);
                
                isGoodCell = ~all(isnan(art1),1) & ~all(isnan(art2),1);
                
                goodCells = find(isGoodCell);
                goodCells = goodCells(randperm(length(goodCells)));
                isGoodCell = false(size(isGoodCell));
                isGoodCell(goodCells(1:matchNum)) = true;
                
                art1 = art1(:,isGoodCell);
                art2 = art2(:,isGoodCell);
                
                xc = corr(art1,art2);
                pearson(si,sj,isGoodCell) = xc(1:length(xc)+1:end);
                pv(si,sj) = corr(art1(:),art2(:));
                
                art1 = rt1(isGood,:);
                art2 = rt2(isGood,:);
                
                art1(isnan(art1)) = 0;
                art2(isnan(art2)) = 0;
                zeroedpv(si,sj) = corr(art1(:),art2(:));
            end
        end     
        out.pearson = pearson;
        out.pv = pv;
        out.zeroedpv = zeroedpv;
        tmp = toc;
        fprintf('  Time:  %0.3fs.',tmp);
    end
end