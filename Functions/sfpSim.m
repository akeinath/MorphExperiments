function sim = sfpSim(SFPs)

    sim = nan([length(SFPs(1,1,:,1)) length(SFPs(1,1,:,1)) length(SFPs(1,1,1,:))]);
    for k = 1:length(SFPs(1,1,1,:))
        for i = 1:length(SFPs(1,1,:,1))
            m1 = SFPs(:,:,i,k);
            if ~all(isnan(m1(:)) | m1(:)==0)
                for j = 1:length(SFPs(1,1,:,1))
                    m2 = SFPs(:,:,j,k);
                    if ~all(isnan(m2(:)) | m2(:)==0)
                        sim(i,j,k) = nansum(m1(:)>0 & m2(:)>0)./ ...
                            nanmin(nansum(m1(:)>0),nansum(m2(:)>0));
                    end
                end
            end
        end
    end
end