
function extract_explosions(varargin)

% Calculates a host of acoustic characteristics of explosions.
%
%
% If called and not simply ran:
% Required Input Arguments:
% 'dfn', string - File name of _ALL file containing 
%        explosion information
% 'dpn', string - Path name of _ALL file containing 
%        explosion information
% 'xw', string - Path name of the XWAV drive containing undecimated
%       acoustic data
% 'tf', string - File name of the transfer function associated with 
%       the deployment being analyzed
% 'tfp', string - Path name of the transfer function
% 'winTE', int - Length of extension of the window used in the Teager
%          energy calculation
% 'tTE', int - Threshold used in Teager energy calculation
%
% If just ran, the program will prompt you for the _ALL file, the XWAV
% drive, and the transfer function file
%
%
% Example: 
%
% extract_explosions('dfn', 'SOCAL30J_DALL.mat', 'dpn',...
%     'J:\Verified Detections\Site J\ALL detectionfiles', ...
%     'xw', 'H:\', 'tf', '483_081205_invSensit.tf', 'tfp',...
%     'J:\Transfer Functions\400_series\483_081205', 'winTE', 3000, 'tTE', 42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Input Arguments and Variable Assignment %%
N = 20000; % make this equal to the decimated sampling rate if you want 1 Hz bins
decFact = 10;
pflag = 0; % plot flag, set to 1 for plots, 0 for no plots
teagerhilb = 0; %include Teager/Hilbert calculations: 1; exclude = 0
headerFileName = [];
percHil = .05;
histFlag = 1; % do you want to make/save histograms?

n = 1;

while n<=length(varargin)
    switch varargin{n};
        case 'dfn' %combined file name and path
            detFileName = varargin{n+1};
            n=n+2;
        case 'dpn'
            detPathName = varargin{n+1};
            n=n+2;
        case 'xw' %xwave path name
            xwavDriveName = varargin{n+1};
            n=n+2;
        case 'tf' %transfer function file name and path
            tfName = varargin{n+1};
            n=n+2;
        case 'tfp'
            tfPathName = varargin{n+1};
            n=n+2;
        case 'winTE' %Teager energy window and threshold
            winTE = varargin{n+1};
            n=n+2;
        case 'tTE'
            thresTE = varargin{n+1};
            n=n+2;
        case 'smooth' %smooth window
            smoothWin = varargin{n+1};
            n=n+2;
        otherwise
            error('Bad argument');
    end
end

% set verified combined explosion detections path name
if ~exist ('detFileName','var')
    [detFileName,detPathName,~] = uigetfile('*.mat','Select verified detection file',...
    'I:\Promotion Anna Stand 20.05.16\Projekt\Verified Explosions Matlab-files\Site H\ALL detections');
end

% choose xwave drive
if ~exist('xwavDriveName','var')
    xwavDriveName = uigetdir('','Select xwave drive');
end

% choose transfer function
if ~exist('tfName', 'var')
    [tfName,tfPathName,~] = uigetfile('*.tf','Select transfer function file','D:\Tethys database\HARP_TF\700_series');
end

if ~exist('winTE', 'var')
   winTE = 300; % Teager Energy smooth used for duration calculation
end

if ~exist('thresTE', 'var')
    thresTE = 28; % Threshold for Teager Engery duration calculation in dB (no transfer function)
end

% if smooth window wasn't passed in, make it 3000
if ~exist ('smoothWin', 'var')
    smoothWin = 300;
end

%% Load data %%
load(fullfile(detPathName,detFileName),'bt');

% remove false positive explosion detections
goodRows = bt(:,3) == 1;
btPruned = bt(goodRows,:);

headerFileName = strrep(detFileName, '.mat', '_headerInfo.mat');

if ~exist(fullfile(detPathName,headerFileName),'file')
   
    % Make a list of xwavs to compare times to
    dirList = dir(fullfile(xwavDriveName,'*disk*'));
    xwavPathAll = [];
    for iD = 1:length(dirList)
        xwavNameList = dir(fullfile(xwavDriveName,dirList(iD).name,'*x.wav'));
        xwavNameMat = vertcat(xwavNameList(:).name);
        xwavPath = fullfile(xwavDriveName,dirList(iD).name);
        xwavPathMat = repmat([xwavPath,'\'],size(xwavNameMat,1),1);
        xwavFullfile = cellstr([xwavPathMat,xwavNameMat]);
        xwavPathAll = [xwavPathAll;xwavFullfile];
    end
    
    nXwav =  size(xwavPathAll,1);
    for iF = 1:nXwav
        [rawStart,rawDur,fs] = readxwavhd(xwavPathAll{iF});
        fileStart(iF,1) = datenum(rawStart(1,:));
        fileEnd(iF,1) = datenum(rawStart(end,:))+ (rawDur(end)/(60*60*24));
        
        if rem(iF,100)==0
            fprintf('Done with file %d of %d \n',iF,nXwav)
        end
    end
    cd(detPathName);
    save(fullfile(detPathName, headerFileName),'fileStart','fileEnd','fs','xwavPathAll')
else
    load(fullfile(detPathName,headerFileName))
end

%% Declaring variables/pulling data from file %%

% get binsize using sampling rate and decmiation factor
fsDF = fs/decFact;
binSize = fsDF/N;

Fc1 = 50;   % High Pass cut off
FO = 10;     % Order
% [B,A] = butter(FO/2,Fc1/(fs/2),'high');

% highpass filter
[b,a] = ellip(4,0.1,40,Fc1*2/fsDF, 'high');


maxFofInt = 10000; % TODO
fHz = Fc1:binSize:(maxFofInt);
[fHzInterp, tfVals] = dtf_map(fullfile(tfPathName,tfName),fHz);
fkHz = fHzInterp/1000;
wExp = zeros(1,N);
spExpMat = [];



binsOfInt = (Fc1/binSize):binSize:(maxFofInt)/binSize;

nExp = size(btPruned,1); % number of explosions to process
peakFreq = zeros(nExp,1); % preallocate vector space
ppSignal = zeros(nExp,1);
spExpMatTf = zeros(nExp,length(binsOfInt));
dB3band = zeros(nExp,1);
dB10band = zeros(nExp,1);
exStart = zeros(nExp, 1); %raw start times of each explosion
IEI = zeros(nExp - 1, 1); %inter explosion intervals

startBuff = (0.2*fs) + 10000; % extra samples are to remove after filtering
endBuff = 0.5*fs;

%% Iterate through good detections   %%
for iE = 1:nExp
  
    % pull data from xwav drive
    % find matching xwav by time
    exStart(iE, 1) = btPruned(iE,4);
    exEnd = btPruned(iE, 5);
    
    % inter-explosion interval
    if iE >= 2   
        IEI(iE-1, 1) = (exStart(iE, 1)-exStart(iE-1, 1))*24*60*60;
    end
    

    tDiff = fileStart - exStart(iE, 1);
    fIdx = find(tDiff<=0,1,'last');
    
    matchFile = xwavPathAll{fIdx};
    
    % read heaader details
    PARAMS = read_header(matchFile);
    
    % identify raw file that starts right before detection
    rawDiff = PARAMS.dnumStart - exStart(iE, 1);
    rawIdx = find(rawDiff<=0,1,'last');
    
    % figure out what byte that is
    rawBytes = PARAMS.byte_loc(rawIdx);
    
  %% Explosion Timing etc. (undecimated) %%
  
    % calculate the offset between detection and raw file
    skip = -rawDiff(rawIdx)*60*60*24*PARAMS.sample_rate(1)*PARAMS.nch*2;
    
    % get detection location by addtion raw start and offset
    % to get the position in the file from which you must start
    % readin
    rawExpLoc = floor(rawBytes + skip - startBuff*2);
    rawExpLocNoBuff = floor(rawBytes + skip);
    
    % make sure that location is a multiple of 2
    if rem(rawExpLoc,2) > 0
        rawExpLoc = rawExpLoc+1;
    end
    
    % to ensure that we don't begin further than the beginning of the raw
    % file
    [rawExpLoc,sFlag]= max([rawExpLoc,PARAMS.byte_loc(1)]);
    
    if sFlag == 1 % if we were able to back up by the buffer amount,
        % calculate true start by subtracting buffer from bt
        buffOffset = startBuff - 10000;
    else % if we hit the file start, our new detection start must be the earliest data time in the file.
        buffOffset = rawExpLocNoBuff - PARAMS.byte_loc(1);
    end
    
    % calculate detection duration (how much of the file to read from the
    % starting point based on btPruned start/end times
    exDur = exEnd - exStart(iE, 1);
    exDurSampNoBuf = exDur*60*60*24*PARAMS.sample_rate(1);
    exDurSamp = exDur*60*60*24*PARAMS.sample_rate(1)+ (2*endBuff) + buffOffset;
    exDurSamp = min(exDurSamp,PARAMS.byte_loc(end)+PARAMS.byte_length(end));
    
    % compute end byte
    % rawExpLocEnd = rawExpLoc + exDurSamp;
    
    % extract timeseries
    fid = fopen(matchFile,'r');
    fseek(fid,rawExpLoc,'bof');%
    dtype = 'int16';
    % read from the start location for the length of exDurSamp
    DATA = fread(fid,[PARAMS.nch,floor(exDurSamp)],dtype)';
    fclose(fid);

    
    %% Decimated Land from here forward!

    % decimate the data
    decData = decimate(DATA, decFact);
    buffOffset = buffOffset/10;
    exDurSampNoBuf = buffOffset + (exDurSampNoBuf/10);
    exTS = filter(b, a, decData);
    exTS = exTS(1001:end); % exclude filter artifact so it doesn't interfere
    % with any subsequent calculations
     
     
  %% Duration Calculation
 
     dcOffset = mean(exTS);
     exTS_abs = abs(exTS-dcOffset); %take the absolute value of raw energies
     exTS_smooth = fastsmooth(exTS_abs, smoothWin); %smooth

     exTS_dB = 20*log10(exTS_smooth); %convert to dB
     exTS_dB_noInf = exTS_dB(~isinf(exTS_dB));
     
     % TODO break here to check if any -infs are being added to the time
     % series
     if length(exTS_dB_noInf) == length(exTS_dB)
         1;
     end
     
     exTS_dB_75 = prctile(exTS_dB_noInf,75); %distribution of RLs

     %exTS_dB_sorted = sort(exTS_dB_noInf);
   
     avg = mean(exTS_dB_noInf(exTS_dB_noInf<=exTS_dB_75)); %take avg over the first 75% of sorted RLs

     [C1,I1]=max(exTS_dB(buffOffset:floor(exDurSampNoBuf))); %find peak intensity within
     % detector's window
     
     I1 = I1 + buffOffset; % for graphing purposes
     
     [C, I] = max(exTS_dB); % find peak intensity over window

     half1 = exTS_dB(1:I); % split the time series into halves around the peak energy
     half2 = exTS_dB(I+1:end);

     stExp = find(half1 <= avg,1,'last')+1; %start and end points of explosion
     endExp = find(half2 <= avg, 1, 'first')-1+length(half1); 

     DurAvg(iE, 1) = (endExp-stExp)/fsDF; % in seconds
     
     ends(iE, 1) = endExp;
     starts(iE, 1) = stExp;
     
     exTSMat{iE,1} = exTS;  
     
     %% Signal:Noise Ratio
 
     %exTS_abs = abs(exTS-dcOffset); %take the absolute value of raw energies  
     
     %exTS_abs_noInf = exTS_abs(~isinf(exTS_abs));
     
     %exTS_dB_75 = prctile(exTS_dB_noInf,75); %distribution of RLs
     %exTS_dB_10 = prctile(exTS_dB_noInf, 10);
     
     %noise = exTS_dB_noInf(exTS_dB_noInf <= exTS_dB_75); %estimate the background noise
     %noise = noise(noise > exTS_dB_10);
     
%      noise = [exTS_dB(1:stExp+1); exTS_dB(endExp+1:end)];
%      noise = noise(~isinf(noise));
%      signal = max(exTS_dB); 
%      
%      SNR(iE, 1) = signal - mean(noise); %compute the SNR
%      RL(iE, 1) = max(exTS_dB_noInf); %compute the RL (not in dB)

     noise = [exTS_dB(1:stExp+1); exTS_dB(endExp+1:end)];
     noise = noise(~isinf(noise));
     signal = max(exTS_dB); 
     
     SNR(iE, 1) = signal - mean(noise); %compute the SNR
     RL(iE, 1) = max(exTS_dB_noInf); %compute the RL (not in dB)
     
  %% Rise Time
  
%      sealThresh = 77; %threshold of hearing for a seal underwater at 800 Hz
%      %Reichmuth, C., Holt, M.M., Mulsow, J. et al. J Comp Physiol A (2013) 199: 491. doi:10.1007/s00359-013-0813-y     
%      
%      rtF1 = 700;  %800 Hz band limits
%      rtF2 = 900;
%      
%      [p1,p2] = ellip(4,0.1,40,[rtF1 rtF2]*2/fs); 
%      exTS_rt = filter(p1, p2, exTS); %filter the signal to a 100Hz band around 800 Hz
     
     %get the SPL of each sample in the time series
     %SPL_rt(:,1) = 20*log10(exTS_rt(:).^2)+ tfVals(800);
%      
%      if SPL_rt(stExp) >= sealThresh
%          stRise = stExp;
%      else
%          stRise = find(SPL_rt(stExp:endExp) >= sealThresh);
%      end
%      
%      %TODO
%      endRise = find(SPL_rt(stRise:endExp) >= sealThresh + 90);
%      
%      riseTime(iE, 1) = (endRise-stRise)/fsDF; 
     
     
  %% Compute parameters
    
    % explosion portion of exTS
    expTS = exTS(stExp:endExp);
    
    % spectrum
    winLen = length(expTS);
    sigWin = hann(winLen);

    wExp(1:winLen) = sigWin'.*expTS';
    spExp = 20*log10(abs(fft(wExp,N)));
    sub = 10*log10(PARAMS.sample_rate(1)/N);
    spExpSub = spExp-sub;
                
    %reduce data to first half of spectra
    spExpMat = spExpSub(1,1:N/2);  
    spExpMatTf(iE,:) = spExpMat(binsOfInt)+tfVals;
   
  %% Peak Frequency, 3dB and 10dB Bandwith/Frequencies, Center Frequency
  
    % compute peak frequency
    [peakFreqAmp, peakFreqIdx] = max(spExpMatTf(iE,:));
    peakFreq(iE,1) = fkHz(peakFreqIdx);
    
    % peak to peak amplitude
    ppNoTf = 20*log10(max(expTS) + abs(min(expTS)));
    ppSignal(iE,1) = ppNoTf + tfVals(peakFreqIdx);
    
    % compute 3db and 10db bandwidths
    dB3 = peakFreqAmp-3;
    dB10 = peakFreqAmp-10;
    
    % divide the spectrum into halves around the peak frequency
    spectHalf2 = spExpMatTf(iE,peakFreqIdx:end);
    spectHalf1 = spExpMatTf(iE,1:peakFreqIdx);
    
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
    dB3band(iE,1) =  hi3dBFreq-low3dBFreq;
    dB10band(iE,1) =  hi10dBFreq-low10dBFreq;
    
    % calculate center frequencies
    if hi3dBFreq/low3dBFreq >= 1.1
        centerFreq(iE, 1) = sqrt(hi3dBFreq*low3dBFreq); %http://www.learningaboutelectronics.com/Articles/Center-frequency-calculator.php
    else
        centerFreq(iE, 1) = (hi3dBFreq + low3dBFreq)/2;
    end
    
    %% Teager and Hilbert Durations
    % will only run if asked for
    if teagerhilb == 1
        
        % Hilbert Energy
        
        % calculate Hilbert Energy for each sample
        pre_env_y=hilbert(exTS.'); 
        env_y=sqrt((real(pre_env_y)).^2+(imag(pre_env_y)).^2); %Au 1993, S.178, equation 9-4

        cumEnv_y= cumsum(env_y.^2)./sum(env_y.^2);

        % find start and end times of explosion; calculate duration
        [~,startHilb] = min(abs(percHil - cumEnv_y));
        [~,endHilb] = min(abs((1-percHil) - cumEnv_y));
        DurHilb(iE,1) = (endHilb-startHilb)/fsDF;

        %Teager Energy

        TEnergy=zeros(size(exTS));
        [rows, columns]=size(exTS);
        
        % compute Teager Energy for each sample
        TEnergy(2:rows-1)=exTS(2:rows-1).^2-(exTS(1:rows-2).*exTS(3:end));

        avgTEnergy = zeros(length(TEnergy),1);
        
        for a = 1:length(TEnergy)-winTE
            avgTEnergy(a) = mean(abs(TEnergy(a:a+winTE)));
        end
        
        % convert to dB
        avgTEnergy_dB = 10*log10(avgTEnergy);

        aboveTEnergyThres = find(avgTEnergy_dB >= thresTE);

        % compute start and end times of each explosion
        if ~isempty(aboveTEnergyThres)
            above_Thr = aboveTEnergyThres(1);
            below_Thr = aboveTEnergyThres(end);
            startIdx = above_Thr + round(winTE/2);
            endIdx = below_Thr + round(winTE/2);
        else
            startIdx = 1;
            above_Thr = 1;
            endIdx = length(avgTEnergy_dB);
            below_Thr = length(avgTEnergy_dB);
        end
        
        % duration
        DurTE(iE,1)=(endIdx-startIdx)/fsDF;
        startOffsetTE(iE,1) = startIdx; 
    end   
     
  %% Sound Exposure Level and Sound Pressure Level
  
      %calculate SEL WHICH ONE IS RIGHT
      SEL(iE,1) = 10*log10(sum(expTS.^2)/DurAvg(iE,1))+ tfVals(peakFreqIdx);
      %SEL(iE,1) = 10*log10(sum(expTS.^2))+ tfVals(peakFreqIdx);
      
      %calculate SPLrms
      %SPLrms(iE,1) = 20*log10(sqrt(sum(expTS.^2)))+ tfVals(peakFreqIdx);
      
      SPLrms(iE, 1) = 20*log10(sqrt(sum(expTS.^2)/DurAvg(iE, 1))) + tfVals(peakFreqIdx);
      
      %calculate SEL based on SPLrms
      SELrms(iE,1) = SPLrms(iE,1) + 10*log10(DurAvg(iE,1));
      
  %% Various Plots     

     if pflag == 1 % make some plots!
         
        % duration plot
        figure;
        plot(exTS_dB_sorted)
        figure;
        plot(exTS_dB);
        hold on;
        plot(stExp, exTS_dB(stExp), '*g');
        plot(endExp, exTS_dB(endExp), '*r');
        title('exTS dB');
        hold off;
        filename = strcat('53H_', int2str(iE), '.fig');
        saveas(gcf, filename);
        set(gcf,'visible','off');
        
        figure(11);clf
        
        % plot Teager and Hilbert Energy only if they've been calculated
        if teagerhilb == 1
         
             subplot(2,2,1)
             plot(exTS)
             xlabel('Samples')
             ylabel('Counts')
             title(['StartTime: ' datestr(btPruned(iE,4))])

             subplot(2,2,2)
             plot(env_y)
             hold on
             plot(startHilb,env_y(startHilb),'*g')
             plot(endHilb,env_y(endHilb),'*r')
             xlabel('Samples')
             ylabel('Hilbert Env')

             subplot(2,2,3)

             %plot(TEnergy)
             plot(avgTEnergy_dB,'b')
             hold on

             plot(startIdx,avgTEnergy_dB(startIdx),'sg','MarkerFaceColor','g')
             plot(endIdx,avgTEnergy_dB(endIdx),'sr','MarkerFaceColor','r')
             xlabel('Samples')
             ylabel('Teager Energy')
    
             subplot(2,2,4)
        end
        
        figure;
        plot(fHz(binsOfInt)/1000,spExpMatTf(iE,:))
        xlabel('Frequency (kHz)')
        ylabel('dBpp re 1\muPa per Hz')

        figure(12)
        plot(DATA)
        xlabel('Samples')
        ylabel('Counts')
        title('Raw data without filter')
 
        % duration graph
        figure;
        secs = (1:length(exTS_dB))/20000;
        plot(secs(:), exTS_dB(:));
        hold on;
        plot(starts(iE)/20000, exTS_dB(starts(iE)), '*g');
        plot(ends(iE)/20000, exTS_dB(ends(iE)), '*r');
        1;
        pause
        close();
     end
    if rem(iE, 100) == 0 || C > C1
%         figure;
%         secs = (1:length(exTS_dB))/20000;
%         plot(secs(:), exTS_dB(:));
%         hold on;
%         plot(starts(iE)/20000, exTS_dB(starts(iE)), '*g');
%         plot(I/20000, exTS_dB(I), 'ko');
%         plot(ends(iE)/20000, exTS_dB(ends(iE)), '*r');
%         1;
        
        figure;
        plot(exTS_dB);
        hold on;
        title(num2str(SNR(iE)));
        plot(stExp, exTS_dB(stExp), 'g*');
        plot(endExp, exTS_dB(endExp), 'r*');
        1;
        %plot(I1, C1, 'ko');
        %plot(buffOffset, exTS_dB(buffOffset), 'go');
        %plot(floor(exDurSampNoBuf), exTS_dB(floor(exDurSampNoBuf)), 'ro');
    end
end

% make/save histograms of different parameters
if histFlag == 1
    
    % save in the same directory as the _ALL file originated from
    cd(detPathName);
    
    graphs = {DurAvg, SEL, SPLrms, centerFreq, dB10band,...
        dB3band, peakFreq, ppSignal};
    titles = {'DurAvg', 'SEL', 'SPLrms', 'centerFreq', 'dB10band',...
        'dB3band', 'peakFreq', 'ppSignal'};
    binSize = [20, 25, 25, 50, 30, 30, 50, 25];
    
    % make and save each histogram
    for i = 1:length(graphs)
        hist(graphs{i}, binSize(i));
        filename = strcat(detFileName(6:8), '_', titles{i});
        title(titles{i});
        saveas(gcf, filename);
        close;
    end
end
%% Saving important parameters/characteristics of the deployment

modifiedStr = strrep(detFileName, '.mat', '_params.mat');
outFileName = fullfile(detPathName, modifiedStr);
exTimes = btPruned(:,4:5);

% save whatever variables you choose
save(outFileName,'ppSignal','peakFreq','spExpMatTf','fkHz','fs',...
    'tfVals','exTimes','dB3band','dB10band', 'N', 'SEL', ...
    'Fc1','FO','bt','btPruned','IEI', 'DurAvg', 'RL', 'SNR', 'fsDF', ...
    'smoothWin', 'ends', 'starts','exTSMat', 'SPLrms', 'SELrms','centerFreq', '-v7.3')

