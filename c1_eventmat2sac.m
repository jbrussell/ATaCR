clear;

setup_parameter;

% inpath_uncorr = './ENAM/DATA/EVENT_preprocess/';
inpath_uncorr = '/Users/jrussel/RESEARCH/PROJ_YoungPacificORCA/DATA/EVENTS/IRIS_XX_5.5_detrend_18sta/';
ZP21_index = 3;
channel = 'BHZ';


%% Load data

if tf_op == 1
    corrseis_path = sprintf('%s/CORRSEIS/',OUTdir);
elseif tf_op ==2
    corrseis_path = sprintf('%s/CORRSEISAVTF/',OUTdir);
end

stadirs = dir(fullfile(corrseis_path));
for ista = 1:length(stadirs)
    station = stadirs(ista).name;
    inpath_corr = sprintf('%s/%s/',corrseis_path,station);
    if ~isdir(inpath_corr)
        continue
    end
    filenames_corr = dir(fullfile(inpath_corr,['*.mat']));
    disp(station);
    % Loop over event files
    for iev = 1:length(filenames_corr);
        if ~exist(fullfile(inpath_corr,filenames_corr(iev).name))
            continue
        end
        load(fullfile(inpath_corr,filenames_corr(iev).name))
        corrdata = corrseis(ZP21_index).timeseries;
        
        % Load data headers
%         load(fullfile(sprintf('%s/%s/%s_%s_%s.mat',inpath_uncorr,corrected.params.eventid, corrected.params.eventid, corrected.params.network, corrected.params.station)));
        sacin = rdsac(fullfile(sprintf('%s/%s/%s.%s.%s.%s.sac',inpath_uncorr,corrected.params.eventid, corrected.params.eventid, corrected.params.network, corrected.params.station, channel)));
        disp(corrected.params.eventid);
        for itr = 1 %1:length(traces)
%             if ~isempty(find(strcmp(traces(itr).channel,channel)==1))
                
                if tf_op ==1
                    opath = sprintf('%s/CORRSEIS_SAC/%s/',OUTdir,corrected.params.eventid);
                elseif tf_op ==2
                    opath = sprintf('%s/CORRSEISAVTF_SAC/%s/',OUTdir,corrected.params.eventid);
                end
                if ~exist(opath)
                    mkdir(opath);
                end
                
%                 header = traces(itr);
%                 H.DELTA = 1./header.sampleRate;
%                 H.KCMPNM = header.channel;
%                 H.KNETWK = header.network;
%                 H.STLA = header.latitude;
%                 H.STLO = header.longitude;
%                 H.STEL = header.elevation;
%                 H.EVLA = header.station
                H = sacin.HEADER;
                data = corrseis(ZP21_index).timeseries;
%                 evid =  datestr(traces(itr).startTime,'yyyymmddhhMM');
                evid = corrected.params.eventid;
                fullevid = [evid,num2str(H.NZSEC,'%02d'),num2str(H.NZMSEC,'%03d')];
                startTime = datenum(fullevid,'yyyymmddhhMMSSFFF');
                sac_path = fullfile(sprintf('%s/%s.%s.%s.%s.sac',opath,evid, sacin.HEADER.KNETWK, sacin.HEADER.KSTNM, sacin.HEADER.KCMPNM));
                mksac(sac_path,data,startTime,H);
                
%             end
        end       
        
    end
    
end
    