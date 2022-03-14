function mkPolar(vals,binSize,doColor)
    if nargin < 2 || isempty(binSize)
        binSize = 15;
    end
    
    if ~iscell(vals)
        vals = {vals};
    end
    
    if nargin < 3 || isempty(doColor)
        doColor = 1-inferno(length(vals)+1);
        doColor(1,:) = [];
    end
    
    radii = 1:length(vals);
    peak = 0.8;
    for i = 1:length(vals)
        d = vals{i};
        bd = histc(d,[-360:binSize:720-binSize]-binSize./2);
        bd = nansum([bd(1:(360./binSize)) bd((360./binSize)+1:2.*(360./binSize)) ...
            bd(2.*(360./binSize)+1:3.*(360./binSize))],2);
        
        bd = bd./nanmax(bd);
        
        hold on
        x = [0:binSize:360]-binSize./2;
        for j = 1:length(bd)
            a = [radii(i).*cosd(x(j)) radii(i).*cosd(x(j+1)) ...
                (radii(i) + peak.*bd(j)).*cosd(x(j+1)) (radii(i) + peak.*bd(j)).*cosd(x(j))]';
            
            b = [radii(i).*sind(x(j)) radii(i).*sind(x(j+1)) ... 
                (radii(i) + peak.*bd(j)).*sind(x(j+1)) (radii(i) + peak.*bd(j)).*sind(x(j))]';
            
            patch(b,a,doColor(i,:),'facecolor','none','edgecolor',doColor(i,:));
        end
        axis equal
        axis square
        set(gca,'xlim',[-(nanmax(radii)+1) (nanmax(radii)+1)], ...
            'ylim',[-(nanmax(radii)+1) (nanmax(radii)+1)])
        axis off
    end
    plot([0 0],get(gca,'ylim'),'color',[0.5 0.5 0.5],'linestyle','--');
    plot(get(gca,'xlim'),[0 0],'color',[0.5 0.5 0.5],'linestyle','--');
    
end