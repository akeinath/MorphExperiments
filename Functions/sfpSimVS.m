function [sim vssims vscounts] = sfpSimVS(SFPs,val,minComps,root)

    sim = sfpSim(SFPs);
    [vssims vscounts] = cellwiseMatCorr(val,sim);
    
    if nargin > 3 && ~isempty(root)
        figure
        set(gcf,'position',[50 50 600 300])
        subplot(1,2,1)
        imagesc(squarify(nanmedian(sim,3)))
        caxis([0 1])
        colorbar
        axis square
        axis off
        subplot(1,2,2)
        hist(vssims(vscounts>minComps),[-1:0.05:1])
        set(gca,'xlim',[-1 1])
        saveFig(gcf,[root],[{'pdf'} {'tiff'} {'jpeg'}]);
    end
end