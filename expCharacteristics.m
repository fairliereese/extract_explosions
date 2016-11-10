% @params
% 'dpl', struct
%       btPruned, start, end, path, filename
% 'i', int - the iteration that the loop is on, will be appended to the
%   outfile names
% 'xwav', string - the file path to the xwav drive of interest
% @return
% saves the _PARAMS file that has arrays and matrices full of useful 
% acoustic characteristics regarding the explosions from that deployment
function expCharacteristics(dpl, i, xwav)

global nDrives N tfPath toSave 

% if there's only 1 drive, we don't need to append i to the file name TODO
% don't need this anymore ( except we need the headerinfo b/c header info
% is specific to the drive it's on
if nDrives == 1
    headerName = strcat(dpl.filename, '_headerInfo.mat');
%     paramsName = strcat(dpl.filename, '_params.mat');
else
    headerName = strcat(dpl.filename, '_', num2str(i), '_headerInfo.mat');
%     paramsName = strcat(dpl.filename, '_', num2str(i), '_params.mat');
end

% if the headerInfo file doesn't already exist
if ~exist(fullfile(dpl.path, headerName),'file')
   
    % Make a list of xwavs to compare times to
    dirList = dir(fullfile(xwav,'*disk*'));
    xwavPathAll = [];
    for iD = 1:length(dirList)
        
        % check to see if disk file is actually a directory we care about
        subString = strfind(dirList(iD).name, dpl.filename(1:7));
        if ~dirList(iD).isdir || isempty(subString)
            continue
        end
        
        xwavNameList = dir(fullfile(xwav,dirList(iD).name,[dpl.filename(1:7), '*x.wav']));
        xwavNameMat = vertcat(xwavNameList(:).name);
        xwavPath = fullfile(xwav,dirList(iD).name);
        xwavPathMat = repmat([xwavPath,'\'],size(xwavNameMat,1),1);
        xwavFullfile = cellstr([xwavPathMat,xwavNameMat]);
        xwavPathAll = [xwavPathAll;xwavFullfile];
        1;
    end
    
    % get the start/end times of each of the xwav files
    nXwav =  size(xwavPathAll,1);
    for iF = 1:nXwav 
        [rawStart,rawDur,fs, ~] = readxwavhd(xwavPathAll{iF});
        fileStart(iF,1) = datenum(rawStart(1,:));
        fileEnd(iF,1) = datenum(rawStart(end,:))+ (rawDur(end)/(60*60*24));
        1;
        if rem(iF,100)==0
            fprintf('Done with file %d of %d \n',iF,nXwav);
        end
    end
    
    % change directories and save the header file
    cd(dpl.path);
    save(fullfile(dpl.path, headerName),'fileStart','fileEnd','fs','xwavPathAll');
    
% if the header file exists, just load that data instead
else
    load(fullfile(dpl.path,headerName));
end

%% declare variables/pull data from file

btPruned = dpl.btPruned;
btStart = dpl.start;
btEnd = dpl.end;

% get binsize using sampling rate and decmiation factor
fsDF = fs/10;
binSize = 1; %TODO fsDF/N;

cutoff = 50;   % High Pass cut off

% highpass filter
[b,a] = ellip(4,0.1,40,cutoff*2/fsDF, 'high');


maxFofInt = N/20; % sampling frequency, /2, and decimation factor factored in
fHz = cutoff:binSize:maxFofInt;
[fHzInterp, tfVals] = dtf_map(tfPath,fHz);
fkHz = fHzInterp/1000;
wExp = zeros(1,fsDF); %TODO changed to be decimated sample rate
spExpMat = [];



binsOfInt = (cutoff/binSize):binSize:(maxFofInt)/binSize;
nExp = size(btPruned,1); % number of explosions to process
peakFreq = []; % preallocate vector space
ppSignal = [];
spExpMatTf = [];
dB3band = [];
dB10band = [];
exStart = []; %raw start times of each explosion
IEI = []; %inter explosion intervals
RLnoise = []; % received level of the noise
startBuff = (0.2*fs) + 10000; % extra samples are to remove after filtering
endBuff = 0.5*fs;
iEGood = 0; % another index that's important for removing faulty tethys explosions
fivesecs = 5.7870e-05; % five seconds' duration

%% iterate through good detections
for iE = 1:nExp
    
    % check to make sure no weird durations are counted
    if btPruned(iE, btEnd) - btPruned(iE, btStart) > fivesecs
        continue;
    elseif btPruned(iE, btEnd) == 0
        continue;
    else
        iEGood = iEGood + 1;
        btPrunedGood(iEGood, :) = btPruned(iE, :);
    end
    
    % find matching xwav by time
    exStart(iEGood, 1) = btPruned(iE, btStart);
    exEnd = btPruned(iE, btEnd);
    
    % interexplosion interval
    if iEGood >= 2
        IEI(iEGood-1, 1) = (exStart(iEGood, 1) - exStart(iEGood-1, 1)) *24*60*60;
    end
    
    tDiff = fileStart - exStart(iEGood, 1);
    fIdx = find(tDiff<=0, 1, 'last');
    matchFile = xwavPathAll{fIdx};
    
    % read header details 
    [~,~,~,PARAMS] = readxwavhd(matchFile);
    
    % find file that starts before detection
    rawDiff = PARAMS.raw.dnumStart - exStart(iEGood, 1);
    rawIdx = find(rawDiff<=0, 1, 'last');
    rawBytes = PARAMS.xhd.byte_loc(rawIdx);
    
    %% explosion timing etc. (undecimated)
    
    % offset b/w detection and raw file
    skip = -rawDiff(rawIdx)*60*60*24*PARAMS.xhd.sample_rate(1)*PARAMS.xhd.NumChannels*2;
    
    % get detection location 
    rawExpLoc = floor(rawBytes + skip - startBuff*2);
    rawExpLocNoBuff = floor(rawBytes + skip);
    
    % make sure location is a multiple of 2
    if rem(rawExpLoc, 2) > 0
        rawExpLoc = rawExpLoc + 1;
    end
    
    % don't back off beginning of file
    [rawExpLoc, sFlag] = max([rawExpLoc, PARAMS.xhd.byte_loc(1)]);
    
    if sFlag == 1 % if backed up by buffer amound
        buffOffset = startBuff - 10000;
    else % if hit files start
        buffOffset = rawExpLocNoBuff - PARAMS.xhd.byte_loc(1);
    end
    
    % calculate detection duration based on btPruned
    exDur = exEnd - exStart(iEGood, 1);
    exDurSampNoBuf = exDur*60*60*24*PARAMS.xhd.sample_rate(1);
    exDurSamp = exDur*60*60*24*PARAMS.xhd.sample_rate(1)+ (2*endBuff) + buffOffset;
    exDurSamp = min(exDurSamp,PARAMS.xhd.byte_loc(end)+PARAMS.xhd.byte_length(end));
    
    % extract timeseries
    fid = fopen(matchFile, 'r');
    fseek(fid, rawExpLoc, 'bof');
    
    dtype = 'int16';
    
    % read from the start location for calculated duration
    DATA = fread(fid,[PARAMS.xhd.NumChannels,floor(exDurSamp)],dtype)';
    fclose(fid);
    
    %% decimate and format data
    
    decData = decimate(DATA, 10);
    buffOffset = buffOffset/10;
    exDurSampNoBuf = buffOffset + (exDurSampNoBuf/10);
    exTS = filter(b, a, decData);
    exTS = exTS(1001:end); % exclude filter artifact
    
    %% duration calculation
    
    dcOffset = mean(exTS);
    exTS = exTS - dcOffset;
    exTS_abs = abs(exTS); % rectify
    exTS_smooth = fastsmooth(exTS_abs, 300); % smooth
    exTS_dB = 20*log10(exTS_smooth); % convert to dB
    exTS_dB_noInf = exTS_dB(~isinf(exTS_dB));
    
    exTS_dB_75 = prctile(exTS_dB_noInf, 75); % 75th % cutoff
    
    avg = mean(exTS_dB_noInf(exTS_dB_noInf<=exTS_dB_75)); % take avg over first 75% 
    
    [C, I] = max(exTS_dB); % find peak intensity
    
    % split into halves around peak
    half1 = exTS_dB(1:I);
    half2 = exTS_dB(I+1:end);
    
    % find start and end times of exp
    stExp = find(half1 <= avg, 1, 'last') + 1;
    if isempty(stExp)
        stExp = 1;
    end
    
    endExp = find(half2 <= avg, 1, 'first') - 1 + length(half1);
    if isempty(endExp)
        endExp = length(exTS_dB);
    end
    
    % store duration in seconds
    DurAvg(iEGood, 1) = (endExp - stExp)/fsDF;
    
    % store start and end times
    ends(iEGood, 1) = endExp;
    starts(iEGood, 1) = stExp;
    
    exTSMat{iEGood, 1} = exTS;
    
    %% SNR calculation
    
    % noise is everything not within the explosion from the window
    noise = [exTS_dB(1:stExp+1); exTS_dB(endExp+1:end)];
    noise = noise(~isinf(noise));
    signal = exTS_dB(stExp:endExp);
    
    % max, avg SNR and RL in counts
    SNRpeak(iEGood, 1) = max(signal) - mean(noise);
    SNR(iEGood, 1) = mean(signal) - mean(noise);
    RL(iEGood, 1) = max(exTS_dB);
    RLnoise(iEGood, 1) = mean(noise) + tfVals(find(fHzInterp == 500));
    
    %% extra parameters
    expTS = exTS(stExp:endExp);
    
    % spectrum
    winLen = length(expTS);
    sigWin = hann(winLen);

    wExp(1:winLen) = sigWin'.*expTS';
    spExp = 20*log10(abs(fft(wExp,fsDF)));
    sub = 10*log10(PARAMS.xhd.sample_rate(1)/fsDF);
    spExpSub = spExp-sub;
                
    %reduce data to first half of spectra
    spExpMat = spExpSub(1,1:fsDF/2); %TODO changed to decimated samplr8  
    spExpMatTf(iEGood,:) = spExpMat(binsOfInt)+tfVals;
    
%% peak frequency, 3dB and 10dB bandwith/frequencies, center frequency
  
    % compute peak frequency
    [peakFreqAmp, peakFreqIdx] = max(spExpMatTf(iEGood,:));
%     if iEGood == 124
%         1;TODO
%     end
    peakFreq(iEGood,1) = fkHz(peakFreqIdx);
    
    % peak to peak amplitude
    ppNoTf = 20*log10(max(expTS) + abs(min(expTS)));
    ppSignal(iEGood,1) = ppNoTf + tfVals(peakFreqIdx);
    
    % compute 3db and 10db bandwidths
    dB3 = peakFreqAmp-3;
    dB10 = peakFreqAmp-10;
    
    % divide the spectrum into halves around the peak frequency
    spectHalf2 = spExpMatTf(iEGood,peakFreqIdx:end);
    spectHalf1 = spExpMatTf(iEGood,1:peakFreqIdx);
    
    % compute lower cutoffs for sample number and frequency
    if ~isempty(find(spectHalf1<dB3,1,'last')+1)
       low3dBIdx = find(spectHalf1<dB3,1,'last')+1;
    else
       low3dBIdx = 1;
    end
    
    if ~isempty(find(spectHalf1<dB10,1,'last')+1)
       low10dBIdx = find(spectHalf1<dB10,1,'last')+1;
    else
       low10dBIdx = 1;
    end
    
    % compute lower cutoffs for sample number and frequency
    if ~isempty(find(spectHalf2<=dB3, 1,'first')-1)
        hi3dBIdx = find(spectHalf2<=dB3, 1,'first')-1;
    else
        hi3dBIdx = size(spExpMatTf,2);
    end
    
    if ~isempty(find(spectHalf2<=dB10,1,'first')-1)
        hi10dBIdx = find(spectHalf2<=dB10, 1,'first')-1;
    else
        hi10dBIdx = size(spExpMatTf,2); 
    end
    
    % find the frequencies corresponding to the found indices
    hi10dBFreq = fHz(binsOfInt(peakFreqIdx+hi10dBIdx));
    hi3dBFreq = fHz(binsOfInt(peakFreqIdx+hi3dBIdx));
    low10dBFreq = fHz(binsOfInt(low10dBIdx));
    low3dBFreq = fHz(binsOfInt(low3dBIdx)); 
   
    % find the frequency bands
    dB3band(iEGood,1) =  hi3dBFreq-low3dBFreq;
    dB10band(iEGood,1) =  hi10dBFreq-low10dBFreq;
    
    % calculate center frequencies
    if hi3dBFreq/low3dBFreq >= 1.1
        centerFreq(iEGood, 1) = sqrt(hi3dBFreq*low3dBFreq); %http://www.learningaboutelectronics.com/Articles/Center-frequency-calculator.php
    else
        centerFreq(iEGood, 1) = (hi3dBFreq + low3dBFreq)/2;
    end
    
    %% sound exposure level and sound pressure level
    
    % calculate SEL 
    SEL(iEGood,1) = 10*log10(sum(expTS.^2))+ tfVals(peakFreqIdx);

    % calculate SPLrms
    SPLrms(iEGood, 1) = 20*log10(sqrt(sum(expTS.^2)/DurAvg(iEGood, 1))) + tfVals(peakFreqIdx);

    % calculate SEL based on SPLrms
    SELrms(iEGood,1) = SPLrms(iEGood,1) + 10*log10(DurAvg(iEGood,1));
    
end

% %% histograms TODO take these out
% 
% cd(dpl.path);
% graphs = {DurAvg, SEL, SPLrms, centerFreq, dB10band,...
%     dB3band, peakFreq, ppSignal, RLnoise};
% titles = {'DurAvg', 'SEL', 'SPLrms', 'centerFreq', 'dB10band',...
%         'dB3band', 'peakFreq', 'ppSignal', 'RLNoise'};
% binSize = [20, 25, 25, 50, 30, 30, 50, 25, 25];
% 
% % make/save each one
% for i = 1:length(graphs)
%     hist(graphs{i}, binSize(i));
%     title(titles{i});
%     1;
%     close;
% end
% 
% %% PARAMS file save TODO don't need anymore
% save(fullfile(dpl.path, paramsName), 'ppSignal','peakFreq','spExpMatTf','fkHz','fs',...
%     'tfVals','dB3band','dB10band', 'N', 'SEL', 'exTSMat', 'RLnoise', 'btPrunedGood', ...
%     'btPruned','IEI', 'DurAvg', 'RL', 'SNR', 'SNRpeak', 'fsDF', ...
%     'ends', 'starts', 'SPLrms', 'SELrms','centerFreq', '-v7.3');

%% add onto previous toSave variables

% timing variables
toSave.DurAvg = [toSave.DurAvg; DurAvg];
toSave.starts = [toSave.starts; starts];
toSave.ends = [toSave.ends; ends];
toSave.btPruned = [toSave.btPruned; btPruned];
toSave.btPrunedGood = [toSave.btPrunedGood; btPrunedGood];
toSave.IEI = [toSave.IEI; IEI];

% amplitude/inensity variables
toSave.RL = [toSave.RL; RL];
toSave.RLnoise = [toSave.RLnoise; RLnoise];
toSave.SEL = [toSave.SEL; SEL];
toSave.SELrms = [toSave.SELrms; SELrms];
toSave.SNR = [toSave.SNR; SNR];
toSave.SNRpeak = [toSave.SNRpeak; SNRpeak];
toSave.SPLrms = [toSave.SPLrms; SPLrms];
toSave.ppSignal = [toSave.ppSignal; ppSignal];
toSave.tfVals = tfVals;

% frequency variables
toSave.centerFreq = [toSave.centerFreq; centerFreq];
toSave.peakFreq = [toSave.peakFreq; peakFreq];
toSave.dB10band = [toSave.dB10band; dB10band];
toSave.dB3band = [toSave.dB3band; dB3band];

% other
toSave.exTSMat = [toSave.exTSMat; exTSMat];
toSave.spExpMatTf = [toSave.spExpMatTf; spExpMatTf];
toSave.fkHz = fkHz;
toSave.fs = fs;
toSave.fsDF = fsDF;
end




