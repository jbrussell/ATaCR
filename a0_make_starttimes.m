% Make start time files for daily data downloads from list of events
% JBR - 7/13/18

% Input eventfile
evfile = './ENAM/event_starttimes.txt';

% Output dayfile
dayfile = './ENAM/day_starttimes.txt';

% Number of days prior to event for calculation of spectra
Ndays = 4;

% Load event list
evlist = textread(evfile,'%s');

fid = fopen(dayfile,'w');
for iev = 1:length(evlist)
    evnum = datenum(evlist(iev),'yyyymmddHHMM');
    daynums = flip(evnum - [1:Ndays]);
    for iday = 1:length(daynums)
        day = datestr(daynums(iday),'yyyymmddHHMM');
        fprintf(fid,'%s\n',day);
    end    
end
fclose(fid);

