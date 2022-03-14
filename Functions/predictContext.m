function predictContext(vals,root)

    doPerms = perms(1:length(vals(:,1,1)));
%     doPerms = doPerms(randperm(length(doPerms(:,1))),:);
    outP = [root(find(ismember(root,'/'),1,'last')+1:end) '_5params'];
    if exist([outP '.mat'],'file')==2
        load([outP '.mat']);
    else
        null_MAE = nan([length(vals(1,1,:)) length(doPerms(:,1))]);
    end
    fprintf(['\n\t\t\tFitting context pattern for all perms:  '])
    clearLen = 0;
    for o = (doPerms')
        clc
        ind = find(ismember(doPerms,o','rows'));
        if ~isnan(null_MAE(1,ind))
            continue
        end
        str = sprintf([num2str(ind) ' of ' num2str(length(doPerms(:,1)))]);
        fprintf([repmat('\b',[1 clearLen]) str])
        clearLen = length(str);
        
        sVals = squarify(vals(o,o,:));
        sVals(isnan(vals)) = nan;

        [out null_MAE(:,ind)] = help_predCellDOS(sVals,1);
        
        save(outP,'null_MAE');
    end

%     testParams = [1:3];
%     [out MAE] = help_predCellDOL(vals,testParams);
%     [out MAE2] = help_predCellDOS(vals,[1 2]);
% 
%     doPerms = perms(1:length(vals(:,1,1)));
%     if exist('ContextPredictionNull_DOS.mat','file')==2
%         load('ContextPredictionNull_DOS');
%     else
%         null_MAE = nan([length(vals(1,1,:)) length(testParams)+1 ...
%             length(doPerms(:,1))]);
%     end
%     fprintf(['\n\t\t\tFitting context pattern for all perms:  '])
%     clearLen = 0;
%     for o = (doPerms')
%         clc
%         ind = find(ismember(doPerms,o','rows'));
%         if ~isnan(null_MAE(1,1,ind))
%             continue
%         end
%         str = sprintf([num2str(ind) ' of ' num2str(length(doPerms(:,1)))]);
%         fprintf([repmat('\b',[1 clearLen]) str])
%         clearLen = length(str);
%         
%         sVals = squarify(vals(o,o,:));
%         sVals(isnan(vals)) = nan;
%         [out sMAE] = help_predCellDOL(sVals,testParams);
%         null_MAE(:,1:length(testParams),ind) = sMAE;
%         [out sMAE2] = help_predCellDOS(sVals,1);
%         null_MAE(:,end,ind) = sMAE2;
%         
%         save('ContextPredictionNull_DOS','null_MAE');
%     end
% % % 
% % %     testParams = [0:3];
% % %     [out MAE] = help_predCellDOL(vals,testParams);
% % %     [out MAE2] = help_predCellDOS(vals,[1]);
% % %     
% % % % % % %     mk Shuffled values within cell % shuffle rows/cols together randomly
% % %     
% % %     sVals = vals;
% % %     for k = 1:length(vals(1,1,:))
% % %         o = randperm(length(vals(:,1,1)));
% % %         tmp = squarify(vals(o,o,k));
% % %         tmp(isnan(vals(:,:,1))) = nan;
% % %         sVals(:,:,k) = tmp;
% % %     end
% % %     [sout sMAE] = help_predCellDOL(sVals,testParams);
% % %     [sout sMAE2] = help_predCellDOS(sVals,[1]);
% % % %     
% % % %     save('ContextPredictionActual','ae');
% % %     
% % %     ae = [];
% % %     MAE = [MAE MAE2];
% % %     sMAE = [sMAE sMAE2];
% % %     for i = 1:length(MAE(1,:))
% % %         ae = [ae MAE(:,i) sMAE(:,i)];
% % %     end
% % %     
% % %     pairwiseNonparametrics(ae,[root '/MAE_']);
% % %     
% % %     close all
% % %     figure(1)
% % %     set(gcf,'position',[50 50 300 225])
% % %     mkBow(log10(ae))
% % %     ylabel('log_1_0(MAE)')
% % %     xlabel('Order of DOL fit')
% % %     saveFig(gcf,[root '/MAE_log'],[{'tiff'} {'pdf'}])
    
%     %%%%% Surrogate data
% 
%     for surrogateParams = testParams
%         [params error] = fitDriftPattern(vals,surrogateParams);
%         fakeV = nan(size(vals));
%         x = [1:length(vals(:,1,1))];
%         for k = 1:length(vals(1,1,:))
%             fakeV(:,:,k) = DOL(x,params(k,:));
%         end
%         fakeV(isnan(vals)) = nan;
%         fakeV = fakeV+randn(size(fakeV)).*0.05; % add a little noise
% 
%         sfVals = fakeV;
%         for k = 1:length(vals(1,1,:))
%             o = randperm(length(fakeV(:,1,1)));
%             tmp = squarify(fakeV(o,o,k));
%             tmp(isnan(fakeV(:,:,1))) = nan;
%             sfVals(:,:,k) = tmp;
%         end
% 
%         [out MAE] = help_predCellDOL(fakeV,testParams);
%         [sout sMAE] = help_predCellDOL(sfVals,testParams);
%         
%         ae = [];
%         for i = 1:length(testParams)
%             ae = [ae MAE(:,i) sMAE(:,i)];
%         end
%     
%         close all
%         figure(1)
%         set(gcf,'position',[50 50 300 225])
%         mkBow(log10(ae))
%         set(gca,'ylim',[-2 0.5])
%         ylabel('log_1_0(MAE)')
%         xlabel('Order of DOL fit')
%         saveFig(gcf,[root '/SurrogateData/MAE_log_SurrogateParams_' ...
%             num2str(surrogateParams)],[{'tiff'} {'pdf'}])
%     end
end











































