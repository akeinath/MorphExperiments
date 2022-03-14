function mkBatchPreload(paths,forceReload)
    
    clc
    close all
    drawnow
    
    if nargin < 2
        forceReload = true;
    end

    warning off all
     %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    labels = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    
    pause_thresh = -2;
    envSize = [17 17];
    
    allPropAligned = [];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        
        slashInds = find(ismember(upiece{mi},'/'));
        top = ['MatlabData/AnalysisPreloads/' upiece{mi}(slashInds(1)+1:end)];
        
        if exist([top '.mat'],'file')~=2 || forceReload
            s = load(paths{isM(1)});
            doAl = help_getAlignmentID(s.alignment,length(isM),paths(isM));
            alignMap = s.alignment(doAl).alignmentMap;
            am = repmat({[]},[1 length(sessions)]);
            uP = repmat({[]},[1 length(sessions)]);
            aGT = repmat({[]},[1 length(sessions)]);
            atm = repmat({[]},[1 length(sessions)]);
            amfr = repmat({[]},[1 length(sessions)]);
            isPC = repmat({[]},[1 length(sessions)]);
            SHCs = repmat({[]},[1 length(sessions)]);
            SICs = repmat({[]},[1 length(sessions)]);
            PFSs = repmat({[]},[1 length(sessions)]);
            PFRs = repmat({[]},[1 length(sessions)]);
            MFRs = repmat({[]},[1 length(sessions)]);
            SNRs = repmat({[]},[1 length(sessions)]);
            aSamp = nan([envSize length(sessions)]);
            SFPs = s.alignment(doAl).regDetails.spatial_footprints_corrected;
            simVecs = repmat({[]},[1 length(sessions)]);
            envs = [];
            tic

            slashInds = find(ismember(paths{1},'/'));
            root = ['Plots/Summary' paths{1}(slashInds(1):slashInds(2)-1)];
            slashInds = find(ismember(upiece{mi},'/'));
            root = [root '/' upiece{mi}(slashInds(end)+1:end) '/ContinuousAnalyses'];
            close all
            drawnow;
            fprintf(['\t\tPreloading Data... '])
            for si = 1:length(sessions)
                s = load(sessions{si},'processed','exclude');
                slashInds = find(ismember(sessions{si},'/'));
                gT = s.processed.trace;
                if isfield(s.processed,'exclude')
                    gT = gT(s.processed.exclude.SFPs,:);
                end
                
                amfr{si} = nanmean(gT,2);
                v = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
                [m os unm] = mkTraceMaps(s.processed.p,gT,v>=pause_thresh,envSize);
                am{si} = m;
                aGT{si} = gT;
                uP{si} = s.processed.p;
                
                aSamp(:,:,si) = os./30;
                
                %%%% compute spatial info
                unm = reshape(unm,[numel(unm(:,:,1)) length(unm(1,1,:))]);
                normS = os(:)./nansum(os(:));
                normM = bsxfun(@rdivide,30.*unm,30.*nanmean(gT(:,v>=pause_thresh),2)');
                SICs{si} = nansum(bsxfun(@times,normS,30.*unm.*log2(normM)),1)';

                tm = m./repmat(nanmax(nanmax(m,[],1),[],2),size(m(:,:,1)));
                PFSs{si} = permute(nansum(nansum(tm>0.75,1),2),[3 1 2]);
                PFRs{si} = permute(nanmax(nanmax(m,[],1),[],2),[3 1 2]).*30;
                isPC{si} = s.processed.splithalf.wholemap_unmatched.p;
                MFRs{si} = nanmean(gT(:,v>=pause_thresh),2);
                if isfield(s.processed,'exclude')
                    isPC{si} = isPC{si}(s.processed.exclude.SFPs,:);
                    SHCs{si} = s.processed.splithalf.wholemap_unmatched.p(s.processed.exclude.SFPs,:);                
                    SNRs{si} = s.processed.snr.whole(s.processed.exclude.SFPs,:);
                end
                envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
            end  
            
            uGT = repmat({[]},[1 length(sessions)]);
            minVecTimes = nanmin(cellfun(@size,simVecs,repmat({2},size(simVecs))));
            umfr = nan([length(alignMap{1}(:,1)) length(sessions)]);
            upfr = nan([length(alignMap{1}(:,1)) length(sessions)]);
            upfs = nan([length(alignMap{1}(:,1)) length(sessions)]);
            usic = nan([length(alignMap{1}(:,1)) length(sessions)]);
            ushc = nan([length(alignMap{1}(:,1)) length(sessions)]);
            usnr = nan([length(alignMap{1}(:,1)) length(sessions)]);
            uIsPC = nan([length(alignMap{1}(:,1)) length(sessions)]);
            um = nan([envSize length(alignMap{1}(:,1)) length(sessions)]);
            uSFPs = nan([length(SFPs{1}(1,:,1)) length(SFPs{1}(1,1,:)) ...
                length(alignMap{1}(:,1)) length(sessions)]);
            for si = 1:length(sessions)
                umfr(alignMap{1}(:,si)~=0,si) = MFRs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                upfr(alignMap{1}(:,si)~=0,si) = PFRs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                upfs(alignMap{1}(:,si)~=0,si) = PFSs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                usic(alignMap{1}(:,si)~=0,si) = SICs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                ushc(alignMap{1}(:,si)~=0,si) = SHCs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                uIsPC(alignMap{1}(:,si)~=0,si) = isPC{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                um(:,:,alignMap{1}(:,si)~=0,si) = am{si}(:,:,alignMap{1}(alignMap{1}(:,si)~=0,si));
                usnr(alignMap{1}(:,si)~=0,si) = SNRs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si));
                uSFPs(:,:,alignMap{1}(:,si)~=0,si) = permute(SFPs{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),:,:),[2 3 1]);
                
                tmp = aGT{si};
                uGT{si} = nan([length(alignMap{1}(:,1)) length(tmp(1,:))]);
                uGT{si}(alignMap{1}(:,si)~=0,:) = aGT{si}(alignMap{1}(alignMap{1}(:,si)~=0,si),:);
            end
            SFPs = compCropSFPs(uSFPs,[31 31]);
            
            clear aGT
            
            checkP(top);
            save(top,'um','umfr','upfr','ushc','envs', ...
                'uIsPC','SFPs','aSamp','uGT','uP','-v7.3');
            
            isGood = ~isnan(umfr);
            propAligned = nan(length(sessions));
            for si = 1:length(sessions)
                for sj = 1:length(sessions)
                    propAligned(si,sj) = nansum(isGood(:,si)&isGood(:,sj)) ./ ...
                        nanmin(nansum(isGood(:,si)),nansum(isGood(:,sj)));
                end
            end
            allPropAligned = [allPropAligned; mat2lag(propAligned)];
            
            clear um uSFPs
        else
            load(top);
        end
    end
end










































