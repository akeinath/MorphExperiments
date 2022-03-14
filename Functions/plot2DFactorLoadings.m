function plot2DFactorLoadings(attributableR2,shuffle_attributableR2)

    boxExtent = 0.25;

    figure(1)
    set(gcf,'position',[50 50 600 250])
    subplot(1,2,1)
    doInclude = ~isnan(attributableR2(:,1));
    scatter(attributableR2(doInclude,1),attributableR2(doInclude,2), ...
        10,hulkcm(attributableR2(doInclude,:),[0 1]).^(1./2));
    hold on
    plot([0 boxExtent],[boxExtent boxExtent],'linestyle','-','color','k');
    plot([boxExtent boxExtent],[0 boxExtent],'linestyle','-','color','k');
    set(gca,'xlim',[0 1],'ylim',[0 1])
    ylabel('Context-attributed r2')
    xlabel('Drift-attributed r2')
    subplot(1,2,2)
    scatter(shuffle_attributableR2(doInclude,1),shuffle_attributableR2(doInclude,2), ...
        10,hulkcm(shuffle_attributableR2(doInclude,:),[0 1]).^(1./2));
    hold on
    plot([0 boxExtent],[boxExtent boxExtent],'linestyle','-','color','k');
    plot([boxExtent boxExtent],[0 boxExtent],'linestyle','-','color','k');
    set(gca,'xlim',[0 1],'ylim',[0 1])
    ylabel('Context-attributed r2')
    xlabel('Drift-attributed r2')
end