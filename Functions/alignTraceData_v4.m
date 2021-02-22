function s = alignTraceData_v4(doPath,s)

    fn = [doPath '/BehavCam_0/timestamps.csv'];    
    fid = fopen(fn);
    dataArray = textscan(fid, '%f%f%f', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
    fclose(fid); 
    bts = [dataArray{1} dataArray{2}];
    
    
    fn = [doPath '/V4c/timestamps.csv'];    
    fid = fopen(fn);
    dataArray = textscan(fid, '%f%f%f', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);
    fclose(fid); 
    mts = [dataArray{1} dataArray{2}];
    
    oldP = s.pos.p;
    oldP = bsxfun(@minus,oldP,nanmin(oldP')');
    oldP = bsxfun(@times,oldP,nanmax(s.environment.size./nanmax(oldP')));
    
    
    unP = s.pos.uninterp;
    unP = bsxfun(@minus,unP,nanmin(unP')');
    unP = bsxfun(@times,unP,nanmax(s.environment.size./nanmax(unP')));
    [unP gaps] = interpNaNs(unP');
    unP = unP';
    unP(:,gaps>=15) = nan;
    
    posFrames = bts;
    traceFrames = mts;
    
    posFrames(posFrames(:,1)<1,:) = [];
    traceFrames(traceFrames(:,1)<1,:) = [];
    
    traceFrames(traceFrames(:,2) >= posFrames(end,2),:) = [];
    traceFrames(traceFrames(:,2) <= posFrames(1,2),:) = [];
    
    linX = linterp(posFrames(:,2),oldP(1,posFrames(:,1)),traceFrames(:,2));
    linY = linterp(posFrames(:,2),oldP(2,posFrames(:,1)),traceFrames(:,2));
    
    unX = linterp(posFrames(:,2),unP(1,posFrames(:,1)),traceFrames(:,2));
    unY = linterp(posFrames(:,2),unP(2,posFrames(:,1)),traceFrames(:,2));
    
%     s.processed.p = [linX'; linY'];
    s.processed.p = [unX'; unY'];
    s.processed.trace = s.calcium.trace(traceFrames(:,1),:)';
    s.processed.validTraceFrames = traceFrames;
    s.processed.posFrames = posFrames;
%     s.properties.partitions = [find(diff(posFrames(:,2))>300) length(posFrames(:,2))];
end