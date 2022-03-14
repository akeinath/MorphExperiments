function [batchVecs group cellInd] = help_mapSVMPrep(maps,isGroup,doDays)

%     tm = maps(:,:,:,doDays);
%     tm = reshape(tm,[numel(tm(:,:,1,1)) length(tm(1,1,:,1)) length(doDays)]);
%     tm = permute(tm,[1 3 2]);
%     
%     batchVecs = reshape(tm,[numel(tm(:,:,1)) length(tm(1,1,:))]);
% 
%     blah = repmat(isGroup(doDays)',[numel(maps(:,:,1,1)) 1]);
%     group = blah(:);
%     
%     exclude = all(isnan(batchVecs),2);
%     batchVecs(exclude,:) = [];
%     group(exclude,:) = [];

    tm = maps(:,:,:,doDays);
    tm = reshape(tm,[numel(tm(:,:,1,1)) length(tm(1,1,:,1)) length(doDays)]);
    cellInd = repmat(1:length(tm(1,:,1)),[length(tm(:,1,1)) 1]);
    tm = reshape(tm,[numel(tm(:,:,1,1)) length(doDays)]);
    cellInd = cellInd(:);
    
    batchVecs = tm;

    group = isGroup(doDays)';
end