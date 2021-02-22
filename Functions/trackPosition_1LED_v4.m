function [pos uninterp] = trackPosition_1LED_v4(doPath)
    folder = dir(doPath);     %Looks at current file location

    totalFrames = 0;        %Total number of frames observed (including previous video itterations)

    % find avi and dat files
    ff = dir([doPath]);
    tmp = {ff(:).name};
    subFolder = [];
    findString = 'BehavCam_';
    for i = 1:length(tmp)
        if length(tmp{i})>=length(findString) && ...
                ismember(tmp{i}(1:length(findString)),{findString})
            subFolder = tmp{i};
        end
    end
    if isempty(subFolder)
        pos = [];
        uninterp = [];
        fprintf('NO BEHAVIOR CAM FOLDER LOCATED.')
    end
    
    aviFiles = dir([doPath '/' subFolder '\*.avi']);
    filePrefix = 'behavCam';
    
    clipNum = nan(1,length(aviFiles));
    for i = 1:length(aviFiles)
        clipNum(i) = str2num(aviFiles(i).name(1:end-4));
    end
    [a b] = sort(clipNum);
    aviFiles = aviFiles(b);
    clipNum = clipNum(b);

    %Find total number of frames and organize vidoes by their timestamps 
    for ind_vid = 1:length(aviFiles)                         %Going through every object in the folder
        video = VideoReader([doPath '/' subFolder '/' aviFiles(ind_vid).name]);   %Read video file              
        totalFrames = totalFrames + video.NumberOfFrames;    %Determine how many frames in file         
    end

    frame  = 1;                 %Frame # of current video
    CurrentFrame = frame;
    pos = nan(2, totalFrames);
    maxVal = nan(1, totalFrames);
    for it = 1 : length(aviFiles)                     %ignore if file is not a video
        videoObject = VideoReader([doPath '/' subFolder '/' aviFiles(ind_vid).name]);   %Read video file      
        numberOfFrames = videoObject.NumberOfFrames;    %Determine how many frames in current video file   

        v = VideoReader([doPath '/' subFolder '/' aviFiles(ind_vid).name]);             %Read video file a gain ( hasFrame() function won't work otherwise)

        fprintf(['\n\t\tTracking Frames (Clip ' num2str(it) '):  '])

        tic
        while hasFrame(v) && frame<= numberOfFrames
            thisFrame = read(videoObject, frame);
%             thisFrame(:,:,[2 3]) = [];
            thisFrame = prod(double(thisFrame),3);
            maxVal(CurrentFrame) = nanmax(thisFrame(:));
            if nanmax(thisFrame(:)) < 150
                frame = frame + 1;
                CurrentFrame = CurrentFrame + 1;
                continue
            end
            [x y] = find(thisFrame == nanmax(thisFrame(:)));
            pos(:,CurrentFrame) = [nanmedian(x) nanmedian(y)]';
            frame = frame + 1;
            CurrentFrame = CurrentFrame + 1;
        end
        toc
        fprintf('\b')
        frame = 1;
    end
    maxVal(maxVal < (2.*10.^6)) = nan;
    uninterp = pos;
    pos = interpNaNs(uninterp')';
    
    figure(1)
    set(gcf,'position',[50 50 1.*nanmax(pos')-nanmin(pos')])
    plot(pos(1,:),pos(2,:),'color',[0.2 0.2 0.2])
    axis equal
    axis off
    qi = find(ismember(doPath,'/'));
    doPath(qi) = '_';
    outP = ['Plots/Paths/' doPath(qi(2)+1:qi(3)-1) '/' doPath(qi(2)+1:end)];
    saveFig(gcf,outP,'tiff')
    close all
end