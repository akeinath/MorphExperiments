function weirdGeometries(paths)
    clc
    close all
    drawnow
    
    pause_thresh = -2;
    envSize = [16 16];
    warning off all
     %% Split by animal
    mouse = [];
    cond = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),2,'last');
        mouse = [mouse; {paths{i}(ind(1)+1:ind(2)-1)}];
        tmp = paths{i}(ind(2)+1:end-4);
        cond = [cond; {tmp(find(ismember(tmp,'_'),1,'last')+1:end)}];
    end
    umouse = unique(mouse);
    
    doConds = [{'Sq1'} {'G1'} {'U1'} {'X1'} {'H1'}];
    for mi = 1:length(umouse)
        fprintf(['\n\tMouse:  ' umouse{mi} '\n'])
        isM = ismember(mouse,umouse(mi));
        
        sessions = paths(isM);
        
        s = load(sessions{1});
        doAl = help_getAlignmentID(s.alignment,length(sessions),sessions);
        alignMap = s.alignment(doAl).alignmentMap{1};
        
        mcond = cond(isM);
        order = [];
        for k = 1:length(mcond)
            order = [order; find(ismember(doConds,mcond(k)))];
        end
        
        sessions = sessions(order);
        alignMap = alignMap(:,order);
        
        root = ['Plots/ContinuousAnalysis/' umouse{mi}];
        am = repmat({[]},[1 length(sessions)]);
        amfr = repmat({[]},[1 length(sessions)]);
        envs = [];
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','exclude');

            gT = s.processed.trace;
            if isfield(s.processed,'exclude')
                gT = gT(s.processed.exclude.SFPs,:);
            end
            
            unv = [0 sqrt(nansum(diff(s.processed.p,[],2).^2,1))].*30;
            unv(unv > 50) = nan;
            unv = interpNaNs(unv')';
            v = conv(unv,fspecial('gauss',[1 200],15),'same');
            amfr{si} = nanmean(gT(:,v>=pause_thresh),2);

            [m os unm] = mkTraceMaps(s.processed.p,gT,v>=pause_thresh,envSize);
            am{si} = m;
            
            envs = [envs; {lower(sessions{si}(find(ismember(sessions{si},'_'),1,'last')+1:end-4))}];
        end
        fprintf(['Done.\n'])

        umfr = nan([length(alignMap(:,1)) length(sessions)]);
        um = nan([envSize length(alignMap(:,1)) length(sessions)]);
        for si = 1:length(sessions)
            umfr(alignMap(:,si)~=0,si) = amfr{si}(alignMap(alignMap(:,si)~=0,si));
            um(:,:,alignMap(:,si)~=0,si) = am{si}(:,:,alignMap(alignMap(:,si)~=0,si));
        end
        
        iter = 0;
        strLength = 0;
        fprintf(['\t\tComputing pairwise map comparisons: '])
        for si = 1:length(sessions)
            for sj = si+1:length(sessions)
                
                iter = iter+1;
                fprintf(repmat('\b',[1 strLength]));
                str = sprintf([ num2str(iter) ' of ' num2str(nchoosek(length(sessions),2))]);
                fprintf(str);
                strLength = length(str);
                
                tmp1 = um(:,:,:,si);
                tmp2 = um(:,:,:,sj);
            
                figure(1)
                set(gcf,'position',[50 50 250.*length(sessions) length(sessions).*250])
                subplot(length(sessions),length(sessions),(si-1).*length(sessions)+sj)
                pvWarp(tmp1,tmp2);
                axis square

            end
        end       
        
    end
    
end