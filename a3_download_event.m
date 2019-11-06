% download_event

% downloads the event data files that will be corrected for tilt and
% compliance noise

% H. Janiszewski 
% hjaniszewski@carnegiescience.edu
% updated 11/17

clear all;

javaaddpath('IRIS-WS-2.0.6.jar');

startlist = './ENAM/event_starttimes.txt'; %'NOISETC_SAMPLE2/eventtimes_ENAMtest.txt'; % list of start times for data download, this will be the beginning of the waveform so specify appropriately, signal of interest should not be at edges of time series
datalength = 6001; %7200; % length of time series after each start time in seconds (default 86400, code not thoroughly tested for other values)

download_networks = 'YO'; % list of networks to download
download_stations = textread('./ENAM/stalist_good.txt','%s'); %'B02B'; % list of stations to download (* for all)
% Channel Names
chz_vec = 'LHZ'; %'HHZ,BHZ'; % list of acceptable names for Z component
ch1_vec = 'LH1'; %'HH1,BH1'; % list of acceptable names for H1 component
ch2_vec = 'LH2'; %'HH2,BH2'; % list of acceptable names for H2 component
chp_vec = 'LDH'; %'HDH,BDH'; % list of acceptable names for P component

datacache = './ENAM/DATA/EVENT/'; %'NOISETC_SAMPLE2/DATA/datacache'; % output folder for data

%%%%% end user input parameters %%%%%

if ~exist(datacache,'dir')
    mkdir(datacache)
end

startlist = textread(startlist,'%s');
chanlist = sprintf('%s,%s,%s,%s',chz_vec,ch1_vec,ch2_vec,chp_vec);

for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   disp(sprintf('Start Time: %s',eventid));
   if length(eventid) == 12
        otime = datenum(eventid,'yyyymmddHHMM'); % look at this when I get back
   elseif length(eventid) == 14
       otime = datenum(eventid,'yyyymmddHHMMSS');
       eventid = eventid(1:12); % for naming purposes only, start time will still be saved to the second in traces file
   end
   starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
   endtime = datestr(otime+datalength/3600/24,'yyyy-mm-dd HH:MM:SS');
   
   stations_info = irisFetch.Stations('channel',download_networks,download_stations,'*',chz_vec,'startTime',starttime,'endTime',endtime);
   
   
   for ista =1:length(stations_info)
       error = 0;
       stnm = stations_info(ista).StationCode;
       network = stations_info(ista).NetworkCode;
       if ~exist(fullfile(datacache,eventid),'dir')
           mkdir(fullfile(datacache,eventid));
       end
       sta_filename = fullfile(datacache,eventid,[eventid,'_',network,'_',stnm,'.mat']);
       if exist(sta_filename,'file')
           disp(['Exist: ',sta_filename,', Skip!']);
           continue;
       end
       disp(['Downloading station: ',stnm,' From:',starttime,' To:',endtime]);
		try
            traces = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime,'includePZ');
			save(sta_filename,'traces');
		catch e
            e.message;
            error = 1;
        end
        if error ==1
            try % to try and get around the missing zeros for some pressure components                
                traces = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime);
                save(sta_filename,'traces');
            catch e
                e.message;
                continue;
            end
        end
    end
   
end