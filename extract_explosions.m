function extract_explosions(varargin)

% Calculates a host of acoustic characteristics of explosions.
%
%
% Required Input Arguments:
% 'xwav1', string - Path name of the first XWAV drive containing undecimated
%       acoustic data
% 'fromTeth', boolean - true will make bt indexing variables into 1, 2;
%       false will make bt indexing variables into 4, 5.
% 'tf', string - File name of the transfer function associated with 
%       the deployment being analyzed
% 'tfp', string - Path name of the transfer function
% 'sampler8', int - sample rate of the acoustic data
%
% Optional Input Arguments:
% 'site', string - Site name, needed if from tethys
% 'project', string - Project name, needed if from tethys
% 'dpln', int - deployment number, needed if from tethys
% 'dfn', string - File name of _ALL file containing 
%        explosion information, needed if NOT from tethys
% 'dpn', string - Path name of _ALL file containing 
%        explosion information, needed if NOT from tethys
% 'xwav2', string - Path name of the second XWAV drive containing undecimated
%       acoustic data, needed if deployment was split up between multiple
%       hard drives
% 'xwav3', string - Path name of the third XWAV drive containing undecimated
%       acoustic data, needed if deployment was split up between multiple
%       hard drives
% 'saveLoc', string - Path name of where you want stuff to save


%% arg parsing & initialization
n = 1;
global nDrives N tfPath toSave 
while n<=length(varargin)
    switch varargin{n};
        case 'dfn' %combined file name and path
            detFileName = varargin{n+1};
            n=n+2;
        case 'dpn'
            detPathName = varargin{n+1};
            n=n+2;
        case 'xwav1' %xwave path name
            xwav1 = varargin{n+1};
            n=n+2;
        case 'xwav2' %second xwave path name
            xwav2 = varargin{n+1};
            n=n+2;
        case 'xwav3' %third xwave path name
            xwav3 = varargin{n+1};
            n=n+2; 
        case 'tf' %transfer function file name and path
            tfName = varargin{n+1};
            n=n+2;
        case 'tfp'
            tfPathName = varargin{n+1};
            n=n+2;
        case 'fromTeth' %whether or not the deploymeny needs to be pulled 
            %tethys
            fromTeth = varargin{n+1};
            n=n+2;
        case 'project' % project
            project = varargin{n+1};
            n=n+2;
        case 'site' %site
            site = varargin{n+1};
            n=n+2;
        case 'dpln' %deployment number
            dpln = varargin{n+1};
            n=n+2;
        case 'sampler8' %sample rate of the acoustic data
            N = varargin{n+1};
            n=n+2;
        case 'saveLoc' 
            saveLoc = varargin{n+1};
            n=n+2;
        otherwise
            error('Bad argument');
    end
end

%% required input catch

% choose xwave drive
if ~exist('xwav1','var')
    xwav{1} = uigetdir('','Select xwave drive');
end

% choose transfer function
if ~exist('tfName', 'var')
    [tfName,tfPathName,~] = uigetfile('*.tf','Select transfer function file','D:\Tethys database\HARP_TF\700_series');
end

% sampling rate will default to 200,000 Hz
if ~exist ('N', 'var')
    N = 200000;
end

% default save location
if ~exist('saveLoc', 'var')
    saveLoc = 'J:\Extractor Results\extractor_output';
end

tfPath = fullfile(tfPathName, tfName);

%% initialization based on compatibility of specific deployment

% initialize some cell arrays
btPruned = {};
%xwav = {}; TODO

% if the detection file is provided
if ~fromTeth
   
    % choose detection file catch if not from tethys
    if ~exist ('detFileName','var')
        [detFileName,detPathName,~] = uigetfile('*.mat','Select verified detection file');
    end
    
    % load the bt variable
    load(fullfile(detPathName,detFileName),'bt');
    
    % remove false positive explosion detections
    goodRows = bt(:,3) == 1;
    btPruned{1} = bt(goodRows,:);
    
    % btPruned indexing variables
    btStart = 4;
    btEnd = 5;
    
    % create names for the saved files
    saveFile = detFileName;
    
% if we need to create the detection file from tethys
else
     q = dbInit('Server', 'bandolero.ucsd.edu');
     dbSpeciesFmt('Input', 'Abbrev', 'NOAA.NMFS.v1');
     [btPruned{1}, ~, ~] = dbGetDetections(q, 'Project', project, 'Site', site, ...
     'Deployment', dpln, 'SpeciesID', 'Anthro', 'Call', 'Explosion');
    
     % btPruned indexing variables
     btStart = 1;
     btEnd = 2;
     
     % create names for saved files
     saveFile = strcat(project, num2str(dpln), site, '_btPruned');
end
    
% if the detection file needs to be split up
if exist('xwav2','var')
    
      % split detections b/w 3 drives
      if exist('xwav3', 'var')
          [btPruned{1}, btPruned{2}, btPruned{3}] = splitALL3(btPruned{1}, xwav1, ...
              xwav2, btStart);
          xwav{1} = xwav1;
          xwav{2} = xwav2;
          xwav{3} = xwav3;
          nDrives = 3; 
          
      % split detections b/w 2 drives
      else
          % split 2
          [btPruned{1}, btPruned{2}] = splitALL(btPruned{1}, xwav1, btStart);
          xwav{1} = xwav1;
          xwav{2} = xwav2;
          nDrives = 2;
      end

% only one xwav drive
else
    xwav{1} = xwav1;
    nDrives =  1;
end

%% initialize the toSave struct
    
    toSave = struct;
    
    % timing variables
    toSave.DurAvg = [];
    toSave.starts = [];
    toSave.ends = [];
    toSave.btPruned = [];
    toSave.btPrunedGood = [];
    toSave.IEI = [];
    
    % amplitude/intensity variables
    toSave.RL = [];
    toSave.RLnoise = [];
    toSave.SEL = [];
    toSave.SELrms = [];
    toSave.SNR = [];
    toSave.SNRpeak = [];
    toSave.SPLrms = [];
    toSave.ppSignal = [];
    toSave.tfVals = [];
    
    % frequency variables
    toSave.centerFreq = [];
    toSave.peakFreq = [];
    toSave.dB10band = [];
    toSave.dB3band = [];

    % other
    toSave.exTSMat = [];
    toSave.spExpMatTf = [];
    toSave.fkHz = [];
    toSave.fs = 0;
    toSave.fsDF = 0;

%% loop through each xwav drive

for i = 1:nDrives
    
    % create a deployment struct to call characteristics with
    dpl = struct;
    dpl.btPruned = btPruned{i};
    dpl.start = btStart;
    dpl.end = btEnd;
    dpl.path = saveLoc;
    dpl.filename = saveFile;
    
    
 
    expCharacteristics(dpl, i, xwav{i});  
end

%% graphs

cd(saveLoc);
graphs = {toSave.DurAvg, toSave.SEL, toSave.SPLrms, toSave.centerFreq...
    toSave.dB10band, toSave.dB3band, toSave.peakFreq, toSave.ppSignal...
    toSave.RLnoise};
titles = {'DurAvg', 'SEL', 'SPLrms', 'centerFreq', 'dB10band',...
        'dB3band', 'peakFreq', 'ppSignal', 'RLNoise'};
binSize = [20, 25, 25, 50, 30, 30, 50, 25, 25];

% make/save each one
for i = 1:length(graphs)
    hist(graphs{i}, binSize(i));
    title(titles{i});
    1;
    close;
end

%% save data
saveFile = strcat(saveFile, '_params');    
save(fullfile(saveLoc, saveFile), 'toSave', '-v7.3');

    