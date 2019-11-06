%% Main driver for codes

%% Downloading Data (if don't already have it...)
% Make start times for daily data from event file
a0_make_starttimes

% Download daily data
% a1_download_data
a1_sac2mat_data

% Preprocess daily data
a2_daydata_preprocess

% Download events to correct
% a3_download_event
a3_sac2mat_event

% Preprocess event data
a4_eventdata_preprocess

%% Calculate spectral properties and apply corrections

% Calculate daily spectra from day-data
b1_dailystaspectra

% QC anomalous spectra
b2_cleanstaspectra

% Calculate transfer functions
b3_calctransfunc

% Correct event data
b4_correctevent

%% Convert final event files to SAC
c5_eventmat2sac