function model_specificDrift(driftGroup,inRoot,doPlot)



    clc
    close all
    %%%%%%%%%%%%%%% prep params

%     popSize = [400 25 150].*3;
    popSize = [300 25 150].*1;
    ca1Div = [30 30 30].*3;
    
    invout = 15;
    
%     if ~doShuffle
%         invout = [15 10]; %75 25
%     else
%         invout = [15 150];
%     end
    
    numDays = 32;
    maxActive = 0.25;
    envLabels = [{'sq1'} {'g1'} {'g2'} {'g3'} {'sq2'} {'sq3'}];
    envs = [];
    
    root = ['Plots/Model/' inRoot];

    fnum = [0];
    tmp = getFilePaths(['MatlabData/Modeling_GroupedDrift/' inRoot],'mat');
    for i = length(tmp):-1:1
        name = tmp{i}(find(ismember(tmp{i},'/'),1,'last')+1:end-4);

        fnum = [fnum; str2num(name(find(ismember(name,'_'),1,'last')+1:end))];
    end
    fnum = nanmax(fnum)+1;
    root = [root num2str(fnum)];
    
    fprintf(['\n\t\t\tSim: ' num2str(fnum) '\n'])
    
    for k = 1:6:numDays
        envLabels(3:6) = envLabels(randperm(4)+2);
        envs = [envs envLabels];
    end
    envs = envs(1:numDays);
    
    envSize = [48 48];
    cornerSize = [-1 3 6 10 13 17];
    envLabels = [{'sq1'} {'sq2'} {'sq3'} {'g3'} {'g2'} {'g1'}];
    envMasks = true([envSize length(envLabels)]);
    for i = 1:length(envLabels)
        envMasks(1:1+cornerSize(i),end-cornerSize(i):end,i) = false;
        envMasks(end-cornerSize(i):end,1:1+cornerSize(i),i) = false;
    end
    
    %%%%%%%%%%%% build input ratemaps
    
    
    % attractor rate maps
    grm = [{[]} {[]} {[]}];
    
    am = zeros([envSize popSize(2) 2]);
    fl = [randi(envSize(1),popSize(2),1) randi(envSize(2),popSize(2),1) ...
        randi(envSize(1),popSize(2),1) randi(envSize(2),popSize(2),1)];
    for k = 1:popSize(2)
        am(fl(k,1),fl(k,2),k,1) = 1;
        am(fl(k,3),fl(k,4),k,2) = 1;
    end
    am = imfilter(am,fspecial('gauss',[30 30],5),'same');
    am = maxnorm(am);
    
    cam = nan([size(am(:,:,:,1)) numDays]);
    for si = 1:numDays
        mask = envMasks(:,:,ismember(envLabels,envs{si}));
        
        if envs{si}(1)=='g'
            tmp = am(:,:,:,2);
        else
            tmp = am(:,:,:,1);
        end
        
        tmp(repmat(~mask,[1 1 length(tmp(1,1,:))])) = nan;
        
        cam(:,:,:,si) = tmp;
    end
    grm{2} = cam;
    
    clear am
    clear cam
    
    % time cells
    cam = nan([envSize popSize(1) numDays]);
    fl = [randi(envSize(1),popSize(1),1) randi(envSize(2),popSize(1),1)];
    fl = repmat(fl,[1 1 numDays]);
    
    for doGroup = driftGroup
        gInd = find(ismember(envs',doGroup{1}));
        for si = 2:length(gInd)
            
            daysBetween = gInd(si) - gInd(si-1);
            for i = 1:daysBetween
                prior = fl(:,:,gInd(si-1)) + randi(9,popSize(1),2)-5;

                prior(prior<1) = 1;
                prior(prior(:,1)>envSize(1),1) = envSize(1);
                prior(prior(:,2)>envSize(2),2) = envSize(2);
            end
            fl(:,:,gInd(si)) = prior;
        end
    end
        
    for si = 1:numDays
        mask = envMasks(:,:,ismember(envLabels,envs{si}));
        
        am = zeros([envSize popSize(1) 1]);
        for k = 1:popSize(1)
            am(fl(k,1,si),fl(k,2,si),k) = 1;
        end
        
        am = imfilter(am,fspecial('gauss',[30 30],5),'same');
        am = maxnorm(am);
        am(repmat(~mask,[1 1 length(am(1,1,:))])) = nan;
        
        cam(:,:,:,si) = am;
    end
    
    grm{1} = cam;
    
    % Geo cells
    
    cam = nan([envSize popSize(3) length(envLabels)]);
    fl = [randi(envSize(1),popSize(3),1) randi(envSize(2),popSize(3),1)];
    am = zeros([envSize popSize(3) 1]);
    for k = 1:popSize(3)
        am(fl(k,1),fl(k,2),k) = 1;
    end
%     am = imfilter(am,fspecial('gauss',[30 30],5),'same');
    
    geoAsign = randi(2,popSize(3),1);
    for si = 1:length(envLabels)
        mask = envMasks(:,:,si);
        [a b] = warpMaps(am,mask);
        
        c = a;
        c(:,:,geoAsign==2) = b(:,:,geoAsign==2);
        
        c(isnan(c)) = 0;
        c = imfilter(c,fspecial('gauss',[30 30],5),'same');
        c(repmat(~mask,[1 1 length(c(1,1,:))])) = nan;
        
        c = maxnorm(c);
        
        cam(:,:,:,si) = c;
    end
    
    cam2 = nan([size(am(:,:,:,1)) numDays]);
    for si = 1:numDays
        cam2(:,:,:,si) = cam(:,:,:,ismember(envLabels,envs{si}));
    end
    
    grm{3} = cam2;
    
    %%%%%%%%%%%%%%%%%% Setup connectivity and build CA1 maps
    
    arm = cat(3,grm{1},grm{2},grm{3});
    rarm = reshape(arm,[numel(arm(:,:,1,1)) length(arm(1,1,:,1)) length(arm(1,1,1,:))]);
    
    
    weights = zeros(nansum(popSize),nansum(ca1Div));
    ca1Bounds = [0 cumsum(ca1Div)];
    
    for k = 1:length(weights(1,:))
        
        tw = randperm(length(weights(:,1)));
        tw = tw(1:invout(1));
        
        vals = rand(length(tw),1);
        vals = vals./nansum(vals);

        weights(tw,k) = vals;
    end
    
%     for k = 1:nansum(ca1Div)
%         
%         isGroup = (k>ca1Bounds(1:3) & k <= ca1Bounds(2:4));
%         
%         wp = [{randperm(popSize(1))} {popSize(1)+randperm(popSize(2))} ...
%             {popSize(1)+popSize(2)+randperm(popSize(3))}]; 
%         
%         wp = [wp(isGroup) {cat(2,wp{:})}];
%         wp{2} = wp{2}(randperm(length(wp{2})));
%         tw = [wp{1}(1:invout(1)) wp{2}(1:invout(2))];
%         
%         while length(unique(tw))<length(tw)
%             wp{2} = wp{2}(randperm(length(wp{2})));
%             tw = [wp{1}(1:invout(1)) wp{2}(1:invout(2))];
%         end
%         
%         vals = rand(length(tw),1);
%         vals = vals./nansum(vals);
%         
%         weights(tw,k) = vals;
%     end
    
    iter = 0;
    strLength = 0;
    fprintf('\n\tComputing CA1 maps day: ');
    ca1Maps = nan([envSize nansum(ca1Div) numDays]);
    for si = 1:numDays
        iter = iter+1;
        fprintf(repmat('\b',[1 strLength]));
        str = sprintf(num2str(si));
        fprintf(str);
        strLength = length(str);
        
        tm = rarm(:,:,si);
%         tm = poissrnd(tm.*5); %poissrnd(bsxfun(@times,tm,firing_rate'));
        twm = bsxfun(@times,permute(tm,[2 1]),permute(weights,[1 3 2]));
        tm = permute(nansum(twm,1),[2 3 1]);
        tm(all(isnan(rarm(:,:,si)),2),:) = nan;
        
%         [a b] = sort(tm,2);
%         thresh = a(:,round(length(a(1,:)).*(1-maxActive)));
%         tm(bsxfun(@lt,tm,thresh)) = 0;
        rtm = reshape(tm,[envSize, length(tm(1,:))]);
        
        rtm = maxnorm(rtm);
        rtm = rtm-(1-maxActive);
        rtm(rtm<0) = 0;
%         rtm = rtm==1;
        
%         rtm = minmaxnorm(rtm);
        
        
        isBad = isnan(rtm);
        rtm(isBad) = 0;
        rtm = imfilter(rtm,fspecial('gauss',[30 30],5),'same');
        rtm(isBad) = nan;
        
        ca1Maps(:,:,:,si) = rtm; %rtm;
    end
    
    if doPlot
        doCells = 10;
        for k = [1:doCells ca1Div(1)+[1:doCells] ca1Div(1)+ca1Div(2)+[1:doCells]]
            figure(1)
            set(gcf,'position',[50 50 800 900])
            for i = 1:32
                subplot(7,6,i)
                imagesc(ca1Maps(:,:,k,i));
                axis equal
                axis off
                colormap parula
                alpha(double(~isnan(ca1Maps(:,:,k,i))));
            end
            saveFig(gcf,[root '/ExampleCells_Sim_' num2str(fnum) '/Cell_' num2str(k)],[{'pdf'} {'tiff'}])
        end
    end
    
    clear twm rtm tm rarm arm
    
    % save data
    checkP(['MatlabData/Modeling_GroupedDrift/' inRoot 'simulation_data_' num2str(fnum)]);
    save(['MatlabData/Modeling_GroupedDrift/' inRoot 'simulation_data_' num2str(fnum)]);
end
























