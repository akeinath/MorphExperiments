function [excludedCellSim cell_ids embedding] = c2cSim(sim,minComps)
    cellSim = nan(length(sim(1,1,:)),length(sim(1,1,:)));
    
    fprintf('\n\tComputing cell to cell matrix similarity and embedding... ')
    tic
    
    %fisher transform matrices?
    transMats = reshape(sim,[numel(sim(:,:,1)) length(sim(1,1,:))]);
    total_comps = nan(length(sim(1,1,:)));
    for i = 1:length(sim(1,1,:))
        for j = i+1:length(sim(1,1,:))
            isGood = (~isnan(transMats(:,i))&~isnan(transMats(:,j)));
            total_comps(i,j) = nansum(isGood);
            if total_comps(i,j)>1
%                 cellSim(i,j) = corr(transMats(isGood,i),transMats(isGood,j),'type','pearson');
                cellSim(i,j) = nanmean([transMats(isGood,i)-transMats(isGood,j)].^2);
            end
        end
    end
    cellSim = squarify(cellSim);
    total_comps = squarify(total_comps);
    cellSim(total_comps<minComps) = nan;
    
    excludedCellSim = cellSim;
    cell_ids = 1:length(cellSim(1,:));
    while any(nansum(isnan(excludedCellSim))>1)
        [blah isBad] = nanmax(nansum(isnan(excludedCellSim)));
        excludedCellSim(isBad,:) = [];
        excludedCellSim(:,isBad) = [];
        cell_ids(isBad) = [];
    end

    tmp = excludedCellSim;
    tmp(logical(eye(size(tmp)))) = 0; %1
%     tmp = 1-(tmp);
    
    opt.dims = 2;
    opt.display = false;
    embedding = IsoMap(tmp,'k',6,opt);
    embedding = embedding.coords{1}';
    scatter(embedding(:,1),embedding(:,2))
    tmp = toc;
    fprintf('  %0.3fs.',tmp);
end