function h = mkBar(d,xl,doC,varargin)
    if nargin<2 || isempty(xl)
        xl = 1:length(d(1,:));
    end
    
    if ~iscell(d)
        tmp = repmat({[]},[1 length(d(1,:))]);
        for i = 1:length(d(1,:))
            tmp{i} = d(:,i);
        end
        d = tmp;
    end
    
    groups = length(d(:,1));
    ticks = length(d(1,:));
    set(gca,'xlim',[0.5 ticks+0.5],'xtick',[1:ticks],'xticklabel',xl,'fontname','arial',...
        'fontsize',10,'fontweight','normal')
    hold on
    
    w = 0.8;
    wpg = w./groups;
    
    if nargin < 3 || isempty(doC)

%         doC = [0.9 0.3 0.3; 0.3 0.3 0.9; 0.2 0.2 0.2; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
%             0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7;...
%             [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
%             0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7]];
        
        
        doC = [[0.1:0.1:1]' [0.1:0.1:1]' [1:-0.1:0.1]'];
        
        groupColor = doC;
    end


    
    dip = 0.025; 
    th = [];
    jitters = repmat({[]},[size(d)]);
    for i = 1:groups
        for j = 1:ticks
            jitters{i,j} = randn(size(d{i,j})).*0.035;
        end
    end
    for i = 1:groups
        for j = 1:ticks
            
%             if length(d{i,j})>=1
%                 h(i,j) = plot(j-(w./2)+(wpg.*(i-1))+wpg./2+jitters{i,j},d{i,j},'linestyle','none','color',doC(i,:).*0.5+0.5,...
%                     'marker','o','markersize',5,'markerfacecolor',doC(i,:).*0.5+0.5);
%             end
            
        end
    end
    
    for i = 1:groups
        for j = 1:ticks
            
            sVals = sort(d{i,j});
            sVals(isnan(sVals)) = [];
            
            cutoff = [nanmean(sVals)-nanstd(sVals) nanmean(sVals)+nanstd(sVals)];
            cutoff = [-inf inf];

            lqr = sVals(floor(length(sVals).*0.25));
            uqr = sVals(floor(length(sVals).*0.75));
            mqr = sVals(floor(length(sVals).*0.5));
            
            doC = groupColor(mod((i-1),length(groupColor(:,1)))+1,:);
            
            doX = [(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))];
            center = nanmean(doX);
            extent = nanmean(abs(doX-center)).*0.85;
            
            plot(nanmean(doX).*ones(1,2),[max(cutoff(1),sVals(1)) lqr; min(cutoff(2),sVals(end)) uqr]', ...
                'color',doC,'linewidth',1.5,'linestyle','-');
            h(i,j) = patch(extent.*[-1 1 1 0.5 -0.5 -1]+center,...
                [lqr lqr mqr-dip.*range(get(gca,'ylim')) ...
                mqr mqr mqr-dip.*range(get(gca,'ylim'))], ...
                'w','facecolor',doC.*0.5+0.5,'edgecolor',doC,'linewidth',1.5);
            h(i,j) = patch(extent.*[-1 1 1 0.5 -0.5 -1]+center,...
                [uqr uqr mqr+dip.*range(get(gca,'ylim')) ...
                mqr mqr mqr+dip.*range(get(gca,'ylim'))],...
                'w','facecolor',doC.*0.5+0.5,'edgecolor',doC,'linewidth',1.5); %doC./2+0.5
            
% % %             dir = 1;
% % %             h(i,j) = patch([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))],...
% % %                 [0 0 repmat(nanmean(d{i,j}),[1 2])],doC(ceil(i),:),'edgecolor',doC(ceil(i),:),'linewidth',1.5,'facecolor','none');
% % % %             set(h(i,j),'facealpha',0.5);
% % % % % %             plot([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))],repmat(nanmean(d{i,j}),[1 2]),...
% % % % % %                 'color',edgeColor(i,:)','linewidth',1)
% % % % % %             if length(d{i,j})>1
% % % % % %                 plot(repmat(nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),[1 2]),...
% % % % % %                     [nanmean(d{i,j})-dir.*nanstd(d{i,j})./sqrt(sum(~isnan(d{i,j}))) ...
% % % % % %                     nanmean(d{i,j})+dir.*nanstd(d{i,j})./sqrt(sum(~isnan(d{i,j})))], ...
% % % % % %                     'color',edgeColor(i,:),'linewidth',1);
% % % % % %             end
% % %             
% % %             if ~isempty(varargin) && ~isempty(varargin{1})
% % %                 th(end+1) = text(nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) ...
% % %                     (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))]),...
% % %                     0,[' ' varargin{1}{i}],'horizontalalignment','left',...
% % %                     'verticalalignment','middle','fontname','arial','fontsize',9,'rotation',90,...
% % %                     'fontweight','regular','color','k');
% % %             end
        end
    end
end