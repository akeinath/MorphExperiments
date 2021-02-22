function [r2penalty tr2] = help_Maps2Var(um,envs,envLabels)

    %%%%%%%%%%%%%%%%%% REPLICATE ANALYSIS
    sim = nan([32 32 length(um(1,1,:,1))]);
    iter = 0;
    strLength = 0;
    fprintf(['\n\t\tComputing pairwise map comparisons: '])
    for si = 1:32
        for sj = si+1:32

            iter = iter+1;
            fprintf(repmat('\b',[1 strLength]));
            str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(32,2))]);
            fprintf(str);
            strLength = length(str);

            tmp1 = um(:,:,:,si);
            tmp2 = um(:,:,:,sj);

            sim(si,sj,:) = xcorr3transform(tmp1,tmp2,[0 0 0]);
        end
    end    
    
    sim = sim(1:32,1:32,:);
    envs = envs(1:32);
    
    envLabel = envLabels;
    morphPoint = ones(1,20).*4;
    features = [];
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        for k = 1:length(envI)
            features(envI(k),:) = [k i];
        end
    end
    lagMat = abs(bsxfun(@minus,[1:length(sim(:,1,1))],[1:length(sim(:,1,1))]'));
    envMat = abs(bsxfun(@minus,features(:,2),features(:,2)'));

    vec = [];
    for gi = 1:ceil(length(envs)./6)
        vec((gi-1).*6+1:(gi).*6) = morphPoint(gi);
    end
    attractorMat = abs(bsxfun(@minus,[features(:,2)'>=vec(1:length(features(:,2)))],...
        [features(:,2)'>=vec(1:length(features(:,2)))]'));


    RDMs = cat(3,lagMat,attractorMat,envMat);

    %%%% make rdm example figure

    RDMs = RDMs./repmat(nanmax(nanmax(RDMs,[],1),[],2),[size(RDMs(:,:,1))]);
    RDMs = 1-RDMs;
    
    [r2penalty tr2] = classifyCellRDMswInteract(sim,-RDMs,[]);
end