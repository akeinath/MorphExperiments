function modelDecoding()
    acrossLow = prepAcross;
    acrossHigh = prepAcross;

    paths = getFilePaths('MatlabData/Model/Selectivity/Low','mat');
    for p = paths'
        load(p{1},'envs','ca1Maps');
%         tmp = poissrnd(ca1Maps.*10);
%         isBad = isnan(tmp);
%         tmp(isBad) = 0;
%         tmp = imfilter(tmp,fspecial('gauss',[30 30],5),'same');
%         tmp(isBad) = nan;
        sim = getPairwiseMapSim(ca1Maps,'pv');
        acrossLow.contextPrediction.rawPV.actual = cat(3,acrossLow.contextPrediction.rawPV.actual, ...
            decodeContext(sim.pv,envs',6,6));
    end
    root = ['Plots/Model/Selectivity/Low'];
    plotDecoding(acrossLow.contextPrediction.rawPV.actual, ...
        [root '/Decoding/Accuracy_RawPV']);
    
    paths = getFilePaths('MatlabData/Model/Selectivity/High','mat');
    for p = paths'
        load(p{1},'envs','ca1Maps');
%         tmp = poissrnd(ca1Maps.*10);
%         isBad = isnan(tmp);
%         tmp(isBad) = 0;
%         tmp = imfilter(tmp,fspecial('gauss',[30 30],5),'same');
%         tmp(isBad) = nan;
        sim = getPairwiseMapSim(ca1Maps,'pv');
        acrossHigh.contextPrediction.rawPV.actual = cat(3,acrossHigh.contextPrediction.rawPV.actual, ...
            decodeContext(sim.pv,envs',6,6));
    end
    root = ['Plots/Model/Selectivity/High'];
    plotDecoding(acrossHigh.contextPrediction.rawPV.actual, ...
        [root '/Decoding/Accuracy_RawPV']);
    
    figure()
    plot(nanmean(acrossHigh.contextPrediction.rawPV.actual(:,2,:),3))
    hold on
    plot(nanmean(acrossLow.contextPrediction.rawPV.actual(:,2,:),3))
    set(gca,'ylim',[0 1])
end