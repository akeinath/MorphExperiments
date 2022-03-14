function plot3DFactorLoadings(attributableR2,shuffle_attributableR2)

    boxExtent = 0.25;

    figure(1)
    set(gcf,'position',[50 50 1200 500])
    subplot(1,2,1)
    doInclude = ~isnan(attributableR2(:,1));
    plotTernary(attributableR2(doInclude,:),[], ...
        factor3dcm(attributableR2(doInclude,:)), ...
        [{'Drift'} {'Context Group'} {'Shape'}]);
    subplot(1,2,2)
    plotTernary(shuffle_attributableR2(doInclude,:),[], ...
        factor3dcm(shuffle_attributableR2(doInclude,:)), ...
        [{'Drift'} {'Context Group'} {'Shape'}]);
end