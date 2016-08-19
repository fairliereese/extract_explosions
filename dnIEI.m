function dnIEI(varargin)

% Groups each inter explosion interval as happening at day or at night 

% Required Input:
% 
% 'lat', int - Latitude of site
% 'long', int - Longitude of site (can be given in a negative number)
%
%
% Example:
% dnIEI('lat', 33, 'long', -119)
%

vidx = 1;

% grab input arguments
while vidx <= length(varargin)
    switch varargin{vidx}
        case 'lat'
            lat = varargin{vidx+1};
            vidx = vidx+2;
        case 'long' 
            long = varargin{vidx+1};
            vidx = vidx+2;
            if long < 0
                long = 360 - abs(long);
            end
    end
end
% load times from a _PARAMS extracted explosions file
[FileName,PathName,~] = uigetfile('*.mat','Select extracted explosion (PARAMS) file',...
    'Computer');
load(fullfile(PathName,FileName), 'btPruned', 'IEI');

% add a zero to the beginning of the IEI vector so that the start of 

% compute start and end dates of deployment
startD = floor(btPruned(1, 4));
endD = floor(btPruned(length(btPruned), 4)) + 1;

% grab the nighttime boundaries
q = dbInit('Server', 'bandolero.ucsd.edu');
night = dbDiel(q, lat, long, startD, endD);

% make the night matrix into an array so it can form the edges of the bins
% to sort explosions into the day or night category
nights = [];
j = 1;
for i = 1:length(night)
    nights(1, j) = night(i, 1);
    nights(1, j+1) = night(i, 2);
    j = j+2;
end

% day/night binning
[dnCounts, dnInd] = histc(btPruned(:,4), nights);

% initialize day/night vectors for each interexplosion interval
IEId = []; % when any of the explosions involved is during the day
IEIn = []; % when both explosions are occuring at night

% loop through each subsequent pair of explosions
for i = 2:length(dnInd)
    if rem(dnInd(i-1),2) == 0 || rem(dnInd(i),2) == 0
        IEId = [IEId, IEI(i-1)];
    else
        IEIn = [IEIn, IEI(i-1)];
    end
end
1;
end