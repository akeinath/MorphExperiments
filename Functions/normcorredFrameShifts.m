function shifts = normcorredFrameShifts(meanFrames,allPrepped)

    numROIs = 1;
    shifts = zeros(2,length(meanFrames));
    roiMF = meanFrames;
    for mfi = 1:length(meanFrames)
%         figure(1)
%         set(gcf,'position',[50 50 800 800])
        
        ref = normMeanFrame(meanFrames{mfi},3);
        
%         [x y] = meshgrid([1:length(ref(1,:,1))]-length(ref(1,:,1))./2+0.5, ...
%         [1:length(ref(:,1,1))]-length(ref(:,1,1))./2+0.5);
%         d2c = sqrt(x.^2 + y.^2);
%         mask = d2c < nanmin(size(ref(:,:,1))).*0.5;
%         ref(~mask) = nan;
%         imagesc(-ref) %(filter2(hSmall,meanFrames{mfi}) - filter2(hLarge, meanFrames{mfi}))
%         colormap gray
%         axis equal
%         
%         mask = false(size(ref));
%         for i = 1:numROIs
%             rect = getrect(); 
%             roi = uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
%             mask(roi(3):roi(4),roi(1):roi(2)) = true;
%         end
%         
%         ref(~mask) = nan;
        roiMF{mfi} = ref;
        allPrepped{mfi} = permute(allPrepped{mfi},[2 3 1]);
    end
    
    regMF = repmat({[]},[length(meanFrames) length(meanFrames)]);
    regPrepped = repmat({[]},[length(meanFrames) length(meanFrames)]);
    for mfi = 1:length(meanFrames)
        for mfj = [1:mfi-1 mfi+1:length(meanFrames)]
            
            
            mf = (cat(3,roiMF{mfi},roiMF{mfj}));
            [d1,d2,T] = size(mf); bound = 0;
            options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',2, ...
                'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
                'overlap_pre',[50 50 1],'overlap_post',[50 50 1],'max_shift',120,'upd_template',false);
            
            [M,shifts,template] = normcorre_batch(mf,options_nr,(mf(:,:,1)));
            
            
            SHIFTS = shifts(2,1);
            for k=1:size(allPrepped{mfj},3)-1
                SHIFTS = [SHIFTS; shifts(2,1)];
            end
            regPrepped{mfi,mfj} = apply_shifts(allPrepped{mfj},SHIFTS,options_nr,bound/2,bound/2);
            registered_meanFrame = apply_shifts(mf,shifts,options_nr,bound/2,bound/2);
            regMF{mfi,mfj} = registered_meanFrame(:,:,2);
        end
    end
    
    minOverlap = nanmin(nansum(nansum(~isnan(cat(3,roiMF{:})),1),2))./2;
    xcLims = 50;
    allXC = nan(xcLims.*2+1,xcLims.*2+1,length(roiMF),length(roiMF));
    for i = round(length(roiMF)./2)
        for j = [1:i-1 i+1:length(roiMF)]
            allXC(:,:,i,j) = pvxcorr(roiMF{i},roiMF{j},[xcLims xcLims],minOverlap);
            [x y] = find(allXC(:,:,i,j)==nanmax(nanmax(allXC(:,:,i,j))));
            ts = [x-ceil(length(allXC(:,:,1,1))./2) y-ceil(length(allXC(:,:,1,1))./2)];
            shifts(:,j) = ts;
        end
    end
    
    footprint_corr = nan(length(Files_ms));
    disp('Stage 1.5 - Using N0RMCorre for an initial pairwise alignment.')
    for i=1:length(Files_ms)
        for j=i+1:length(Files_ms)

            fp1 = MeanFrames{i}; fp1(fp1<(mean(fp1(:)) - 3*std(fp1(:)))) = mean(fp1(:)); fp1 = fp1 - min(fp1(:)); %fp1 = fp1/mean(fp1(:));
            fp2 = MeanFrames{j}; fp2(fp2<(mean(fp2(:)) - 3*std(fp2(:)))) = mean(fp2(:)); fp2 = fp2 - min(fp2(:)); %fp2 = fp2/mean(fp2(:));
            mf = cat(3,fp1,fp2);
            [d1,d2,T] = size(mf); bound = 0;
            options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
                'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
                'overlap_pre',50,'overlap_post',50,'max_shift',120);

            [M,shifts,template] = normcorre_batch(mf,options_nr,double(mf(:,:,1)));
            SHIFTS = shifts(2,1);
            for k=1:size(Footprints{j},3)-1
                SHIFTS = [SHIFTS; shifts(2,1)];
            end
            registered_footprint = apply_shifts(Footprints{j},SHIFTS,options_nr,bound/2,bound/2);
            registered_meanFrame = apply_shifts(mf,shifts,options_nr,bound/2,bound/2);

            spatial_footprints{1,1} = Footprints{i};
            spatial_footprints{1,2} = registered_footprint;

            meanFrame_Image{1,1} = registered_meanFrame(:,:,1);
            meanFrame_Image{1,2} = registered_meanFrame(:,:,2);

            n = ['sessions_', sprintf('%02d',i), '_', sprintf('%02d',j), '.mat'];
            save([norm_cor_outputs, '/', n], 'spatial_footprints', 'meanFrame_Image');

            s1 = mean(spatial_footprints{1,1},3); s2 = mean(spatial_footprints{1,2},3);
            c = corrcoef(s1(:), s2(:));
            footprint_corr(i,j) = c(1,2); footprint_corr(j,i) = c(1,2);

            disp([i, j, c(1,2)])

        end
    end
end