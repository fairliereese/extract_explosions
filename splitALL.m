% @params
% 'btPruned', matrix - contains start and end times of each explosion
% 'xwav1', string - path name to first xwav drive that has acoustic data
% 'btStart', int - index to the start time of the explosion
% @return
% 'btPruned1', 'btPruned2', matrices - contain start and end time of each
%   explosion, split up by where they occur within the xwav drives
function [btPruned1, btPruned2] = splitALL(btPruned, xwav1, btStart)

    % Make a list of xwavs to compare times to
    dirList = dir(fullfile(xwav1,'*disk*'));
    xwavPathAll = [];
    for iD = 1:length(dirList)
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
    fileEnd = datenum(rawStart(end,:))+ (rawDur(end)/(60*60*24));

    % temp variable to store all explosions
    temp = btPruned;
    
    % explosions on the first xwav drive
    btPruned1 = temp(temp(:, btStart) <= fileEnd, :);

    % explosions on the second xwav drive
    btPruned2 = temp(temp(:, btStart) > fileEnd, :);

end