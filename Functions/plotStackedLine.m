function h = plotStackedLine(val)
    
    toPlot = [];
    for i = 1:length(val(:,1,1))
%         for j = 1:length(val(1,:,1))
            toPlot{i} = permute(val(i,:,:),[3 2 1]);
%         end
    end
    
    doC = (1-v2rgb(1:length(toPlot)+1)).*0.5 + 0.25;
    doC(1,:) = [];
    
    h = mkLine(toPlot(1,:),[],doC);
end






























