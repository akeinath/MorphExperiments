function out = clustParams(params,cInds,distance,chooseK)

    if nargin < 2 || isempty(cInds)
        cInds = 1:length(params(:,1));
    end

    if nargin < 3 || isempty(distance)
        distance = 'euclidean';
    end
    
    if nargin < 4 || isempty(chooseK)
        chooseK = 12;
    end
    
    close all

% % %     c2c_embedding = mdscale(squareform(pdist(params,distance)),2, ...
% % %         'criterion','stress','start','cmdscale');

% % %     [c2c_embedding v] = cmdscale(squareform(pdist(params,distance)));
% % %     c2c_embedding = c2c_embedding(:,1:2);
    
% % %     opt.dims = 2;
% % %     opt.display = false;
% % %     embedding = IsoMap(squareform(pdist(params,distance)),'k',100,opt);
% % %     c2c_embedding = embedding.coords{1}';
        
    [c2c_embedding b c d] = run_umap(params,'n_neighbors',100, ...
        'n_epochs',500, ...
        'distance',distance,'NSMethod','exhaustive','cluster_method_2D','dbm');
    
    ass = kmedoids(c2c_embedding,chooseK,'replicates',5);    


    kVal = nanmax(ass);
    
    
    clear out
    [a b] = sort(ass);
    out.inds = repmat({[]},[1 kVal]);
    for i = 1:kVal
        out.inds{i} = cInds(ass==i);
    end
    [a orderClusters] = sort(cellfun(@length,out.inds),'descend');
    
    out.inds = out.inds(orderClusters);
    newAss = nan(size(ass));
    for i = 1:kVal
        newAss(ass==orderClusters(i)) = i;
    end
    out.assignment = newAss;
    out.params = params;
    out.embedding = c2c_embedding;
    
    tmp = toc;
    fprintf('  %0.3fs.',tmp);
end