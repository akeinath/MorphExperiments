function sim = getPairwiseMapSim(um,varargin)

    if nargin<2 || isempty(varargin)
        varargin = {'pearson'};
    end

    if ismember({'pearson'},varargin)
        %%% Pearson correlation similarity
        fprintf('\n\tComputing pearsons map similarity...')
        tic
        pearson = nan([length(um(1,1,1,:)) length(um(1,1,1,:)) length(um(1,1,:,1))]);
        for si = 1:length(um(1,1,1,:))
            rt1 = reshape(um(:,:,:,si),[numel(um(:,:,1,1)) length(um(1,1,:,1))]);
            for sj = si+1:length(um(1,1,1,:))
                rt2 = reshape(um(:,:,:,sj),[numel(um(:,:,1,1)) length(um(1,1,:,1))]);
                isGood = ~all(isnan(rt1),2) & ~all(isnan(rt2),2); 
                
                xc = corr(rt1(isGood,:),rt2(isGood,:));
                
                pearson(si,sj,:) = xc(1:length(xc)+1:end);
            end
        end     
        sim.pearson = pearson;
        tmp = toc;
        fprintf('  %0.3fs.',tmp);
    end
    
    if ismember({'pv'},varargin)
        %%% Pearson correlation similarity
        fprintf('\n\tComputing population vector map similarity...')
        tic
        pv = nan([length(um(1,1,1,:)) length(um(1,1,1,:))]);
        for si = 1:length(um(1,1,1,:))
            tmp1 = um(:,:,:,si);
            for sj = si+1:length(um(1,1,1,:))
                tmp2 = um(:,:,:,sj);
                
                rt1 = reshape(tmp1,[numel(tmp1(:,:,1,1)) length(tmp1(1,1,:,1))]);
                rt2 = reshape(tmp2,[numel(tmp2(:,:,1,1)) length(tmp2(1,1,:,1))]);
                isGood = ~all(isnan(rt1),2) & ~all(isnan(rt2),2); 
                
                art1 = rt1(isGood,:);
                art2 = rt2(isGood,:);
                
                isGoodCell = ~all(isnan(art1),1) & ~all(isnan(art2),1);
                
                art1 = art1(:,isGoodCell);
                art2 = art2(:,isGoodCell);
                
                pv(si,sj) = corr(art1(:),art2(:));
                
            end
        end     
        sim.pv = pv;
        tmp = toc;
        fprintf('  %0.3fs.',tmp);
    end
    
    if ismember({'cosine'},varargin)
        fprintf('\n\tComputing cosine similarity...')
        tic
        cossim = nan([length(um(1,1,1,:)) length(um(1,1,1,:)) length(um(1,1,:,1))]);
        for si = 1:length(um(1,1,1,:))
            for sj = si+1:length(um(1,1,1,:))
                tmp1 = um(:,:,:,si);
                tmp2 = um(:,:,:,sj);

                tmp1 = reshape(tmp1,[numel(um(:,:,1,1)) length(um(1,1,:,1))]);
                tmp2 = reshape(tmp2,[numel(um(:,:,1,1)) length(um(1,1,:,1))]);

                isGood = ~all(isnan(tmp1),2) & ~all(isnan(tmp2),2);
                tmp1 = tmp1(isGood,:);
                tmp2 = tmp2(isGood,:);

                isGood = ~all(isnan(tmp1),1) & ~all(isnan(tmp2),1);
                tmp1 = tmp1(:,isGood);
                tmp2 = tmp2(:,isGood);

                dp = pdist2(tmp1',tmp2','cosine');
                
                cossim(si,sj,isGood) = diag(dp);
            end
        end     
        sim.cosine = cossim;
        tmp = toc;
        fprintf('  %0.3fs.',tmp);
    end

    if ismember({'peak'},varargin)
        fprintf('\n\tComputing difference in peak locations...')
        tic
        isPeak = um==repmat(nanmax(nanmax(um,[],1),[],2),size(um(:,:,1,1)));
        ploc = nan(2,length(um(1,1,1,:)),length(um(1,1,:,1)));
        d2p = nan(length(um(1,1,1,:)),length(um(1,1,1,:)),length(um(1,1,:,1)));
        for k = 1:length(um(1,1,:,1))
            for i = 1:length(um(1,1,1,:))
                [x y] = find(isPeak(:,:,k,i));
                if ~isempty(x)
                    ploc(:,i,k) = [nanmedian(x) nanmedian(y)]';
                end
            end
            d2p(:,:,k) = sqrt(bsxfun(@minus,ploc(1,:,k),ploc(1,:,k)').^2 + ...
                bsxfun(@minus,ploc(2,:,k),ploc(2,:,k)').^2);
        end
        doToss = repmat(logical(eye(size(d2p(:,:,1)))),[1 1 length(d2p(1,1,:))]);
        d2p(doToss) = nan;
        sim.d2p = d2p;
        tmp = toc;
        fprintf('  %0.3fs.',tmp);
    end

    if ismember({'raw'},varargin)
    
        fprintf('\n\tComputing map distance similarity...')

        normer = nanmax(nanmax(nanmax(um,[],1),[],2),[],4);
        normer = repmat(normer,[size(um(:,:,1,:))]);
        num = um./normer;
        num = reshape(num,[numel(num(:,:,1,1)) length(num(1,1,:,1)) length(num(1,1,1,:))]);

        tic
        map_distance = nan([length(um(1,1,1,:)) length(um(1,1,1,:)) length(um(1,1,:,1))]);
        for si = 1:length(um(1,1,1,:))
            tmp1 = num(:,:,si);
            for sj = si+1:length(um(1,1,1,:))
                tmp2 = num(:,:,sj);

                map_distance(si,sj,:) = nanmean(abs(tmp1-tmp2),1);
            end
        end     
        sim.map_distance = map_distance;
        tmp = toc;
        fprintf('  %0.3fs.',tmp);
    end
end