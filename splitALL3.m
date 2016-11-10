% @params
% 'btPruned', matrix - contains start and end times of each explosion
% 'xwav1', string - path name to first xwav drive that has acoustic data
% 'btStart', int - index to the start time of the explosion
% 'saveFile', contains name of deployment, used to check to see if accessed
%   disk elements are correct
% @return
% 'btPruned1', 'btPruned2, btPruned3', matrices - contain start and end time of each
%   explosion, split up by where they occur within the xwav drives
function [btPruned1, btPruned2, btPruned3] = splitALL3(btPruned, xwav1, xwav2, btStart, saveFile)

    % Make a list of xwavs from first drive to compare times to
    dirList = dir(fullfile(xwav1,'*disk*'));
    xwavPathAll = [];
    for iD = 1:length(dirList)
        
        % check to see if disk file is actually a directory we care about
        subString = strfind(dirList(iD).name, saveFile(1:7));
        if ~dirList(iD).isdir || isempty(subString)
            continue
        end
        
        xwavNameList = dir(fullfile(xwav1,dirList(iD).name,'*x.wav'));
        xwavNameMat = vertcat(xwavNameList(:).name);
        xwavPath = fullfile(xwav1,dirList(iD).name);
        xwavPathMat = repmat([xwavPath,'\'],size(xwavNameMat,1),1);
        xwavFullfile = cellstr([xwavPathMat,xwavNameMat]);
        xwavPathAll = [xwavPathAll;xwavFullfile];
    end

    iF =  size(xwavPathAll,1);
    [rawStart,rawDur,fs] = readxwavhd(xwavPathAll{iF});
    fileStart = datenum(rawStart(1,:));
    fileEnd1 = datenum(rawStart(end,:))+ (rawDur(end)/(60*60*24));

    % Make a list of xwavs from second drive to compare times to
    dirList = dir(fullfile(xwav2,'*disk*'));
    xwavPathAll = [];
    for iD = 1:length(dirList)
        xwavNameList = dir(fullfile(xwav2,dirList(iD).name,'*x.wav'));
        xwavNameMat = vertcat(xwavNameList(:).name);
        xwavPath = fullfile(xwav2,dirList(iD).name);
        xwavPathMat = repmat([xwavPath,'\'],size(xwavNameMat,1),1);
        xwavFullfile = cellstr([xwavPathMat,xwavNameMat]);
        xwavPathAll = [xwavPathAll;xwavFullfile];
    end

    iF =  size(xwavPathAll,1);
    [rawStart,rawDur,fs] = readxwavhd(xwavPathAll{iF});
    fileStart = datenum(rawStart(1,:));
    fileEnd2 = datenum(rawStart(end,:))+ (rawDur(end)/(60*60*24));

    % first drive's detections
    btPruned1 = btPruned(btPruned(:, btStart) <= fileEnd1, :);

    % second drive's detections
    btPruned2 = btPruned(btPruned(:, btStart) > fileEnd1, :);
    btPruned2 = btPruned2(btPruned2(:, btStart) <= fileEnd2, :);

    % third drive's detections
    btPruned3 = btPruned(btPruned(:, btStart) > fileEnd2, :);
    
end
