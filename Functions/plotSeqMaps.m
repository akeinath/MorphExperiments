function plotSeqMaps(um,envs,envLabel,SFPs,root)
%     for si = 1:5
%         goodS = (si-1).*6+1:(si).*6+2;
%         tum = um(:,:,:,goodS);
%         
%         tenvs = envs(goodS);
%         [a b] = ismember(tenvs,envLabel);
%         [a order] = sort(b,'ascend');
%         
%         tum = tum(:,:,:,order);
%         numMissing = nansum(all(all(isnan(tum),1),2),4);
%         doC = find(numMissing==0);
%         helpPlotMaps(tum(:,:,doC,:),[root '/Maps/8_days/Sequence_' num2str(si)])
%         
%         tSFPs = SFPs(:,:,goodS,doC);
%         helpPlotSFPs(tSFPs,[root '/SFPs/8_days/Sequence_' num2str(si)])
%     end
    
    mostMissing = 1;

    numMissing = nansum(all(all(isnan(um),1),2),4);
    for k = find(numMissing<=mostMissing)'
        plotStackConMaps(um(:,:,k,:),envs,envLabel)
        saveFig(gcf,[root '/Maps/Stacked_6_ByContext/Cell_' num2str(k)],[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
    
    numMissing = nansum(all(all(isnan(um),1),2),4);
    for k = find(numMissing<=mostMissing)'
        plotStackConMaps(um(:,:,k,:),envs,[])
        saveFig(gcf,[root '/Maps/Stacked_6_ByTime/Cell_' num2str(k)],[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
end