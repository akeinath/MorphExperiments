function plotStackConMaps(m,envs,order)

    seq = floor(length(m(1,1,:))./6);
    figure
    set(gcf,'position',[50 50 6.*150 seq.*150])
    for i = 1:seq
        
        doI = (i-1).*6+1:nanmin((i-1).*6+6,length(m(1,1,:)));
        if nargin < 3 || isempty(order)
            sI = doI;
        else
            [a b] = ismember(envs(doI),order);
            [a b] = sort(b);
            sI = doI(b); %doI(b);
        end
        
        tcm = [];
        for j = sI
            if isempty(tcm)
                tcm = m(:,:,j)./nanmax(nanmax(m(:,:,j)));
            else
                tcm = cat(2,tcm,nan(length(tcm(:,1)),1),m(:,:,j)./nanmax(nanmax(m(:,:,j))));
            end
        end
        subplot(seq,1,i)
        imagesc(tcm)
        colormap parula
        alpha(double(~isnan(tcm)))
        axis equal
        axis off
    end
end