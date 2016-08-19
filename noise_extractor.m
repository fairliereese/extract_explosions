% makes a vector based on start and end times of an explosion as calculated
% by extract_explosions.m

% load times from a _PARAMS extracted explosions file
[FileName,PathName,~] = uigetfile('*.mat','Select extracted explosion (PARAMS) file',...
    'Computer');
load(fullfile(PathName,FileName), 'exTSMat', 'starts', 'ends');

% average RLs of the noise
RLnoise = [];

for i = 1:length(exTSMat)
     
     % load the windowed explosion
     exTS = exTSMat{i};
    
     dcOffset = mean(exTS); % make sure that the time series is centered around 0
     exTS_abs = abs(exTS-dcOffset); % take the absolute value of raw energies
     exTS_smooth = fastsmooth(exTS_abs, 3000); % smooth

     exTS_dB = 20*log10(exTS_smooth); % convert to dB
     
     noise = [exTS_dB(1:starts(i)+1); exTS_dB(ends(i)+1:end)]; % use the start and
     % end indicies to calculate the noise
     
     % exclude any values that might throw the avg RL calculation off
     noise = noise(~isinf(noise));
     
     % calculate the average RL of the noise
     RLnoise(i) = mean(noise);
end
1;

     