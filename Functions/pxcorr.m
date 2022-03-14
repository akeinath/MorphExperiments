function am = pxcorr(m1,m2,lims,minSamples)
    %only takes 2d maps, computes corrs in parallel
    if nargin==1 || isempty(m2)
        m2 = m1;
    end

    if any(size(m2)<size(m1))
        tmp = m2;
        m2 = m1;
        m1 = tmp;
    end

    if any(size(m2)~=size(m1))
        biggest = max(size(m2),size(m1));
        m1 = padarray(m1,floor([biggest-size(m1)]./2),nan,'both');
        m2 = padarray(m2,floor([biggest-size(m2)]./2),nan,'both');
        
        biggest = max(size(m2),size(m1));
        m1 = padarray(m1,([biggest-size(m1)]),nan,'pre');
        m2 = padarray(m2,([biggest-size(m2)]),nan,'pre');
    end
    
    if nargin<4
        minSamples = 0;
    end
    
    if nargin<3 || isempty(lims)
        lims = [size(m2(:,:,1))-1 0];
    end
    
    if length(lims)==1
        lims = [lims lims 0];
    end
    
    
    [a b c] = size(m1);
    maxLims = [(size(m2(:,:,1))-1) 180];
    lims(lims>maxLims)=maxLims(lims>maxLims);
    
    % toss if all nans (DOESN'T WORK IF NONCONTINUOUS, just for speed with roi-xcs)
    tossCol = all(isnan(m1)&isnan(m2),1);
    tossRow = all(isnan(m1)&isnan(m2),2);
    m1(tossRow,:) = [];
    m1(:,tossCol) = [];
    m2(tossRow,:) = [];
    m2(:,tossCol) = [];
    
    m1 = padarray(m1,[lims(1) lims(2)],nan);
    m2 = padarray(m2,[lims(1) lims(2)],nan);
    
    xlims = -lims(1):lims(1);
    ylims = -lims(2):lims(2);
    rlims = -lims(3):lims(3);
    
    am = nan(length(xlims),length(ylims),length(rlims));
    for r = 1:length(rlims)
        rm2 = imrotate(m2,rlims(r),'crop');
        sm2 = imrotate(true(size(m2)),rlims(r),'crop');
        rm2(~sm2) = nan;
        for i = 1:length(xlims)
            xm2 = circshift(rm2,[xlims(i) 0]);
            for j = 1:length(ylims)
                tm2 = circshift(xm2,[0 ylims(j)]);
                isGood = ~isnan(m1)&~isnan(tm2);
                if nansum(isGood(:))>=minSamples
                    am(i,j,r) = corr(m1(isGood),tm2(isGood));
                end
            end
        end
    end
end





















