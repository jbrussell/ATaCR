% parameters seeting for calculating OBS tilt and compliance 
clear all;
close all;

% path for matlab codes and functions
addpath ( './function');
PROJ = 'OUT_COMPLIANCE_day';

% location of the continures sacdata for event based
WORKINGdir = '/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz/'; %'/Users/helenj/Cascadia/EARTHQUAKES/YEAR3/COMPLIANCESED/J42C/';

% EVTsacdata_input, SACPZ_input directory 
INSTRUMENTdir =  'NONE'; %'../INSTRUMENT/';
EVTsacdir     = '/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz/'; %'/Users/helenj/Cascadia/EARTHQUAKES/YEAR3/COMPLIANCESED/J42C/';

stalist = textread(['./stalist_good.txt'],'%s\n');

% OUTPUT PATHS
compliancematpath = ['./',PROJ,'/COMPLIANCEMAT/'];

% channel naming
chz='HZ';            %Channel name of Z component
ch1='H1';    %Channel name of H1 component 
ch2='H2';    %Channel name of H2 component
chp='DH'; %'XH';	  %Channel name of DPG component

% make sure you have same sampling rate for all 4 components 
dt= 1; %1/125;    

% day info to calcautle tilt and compliance
nday = 1; %2;     % how many days before event time.
T    = 6000;  % the legnth of each time window, in sec  
nwin = 20;    % numbers of time window in a day
% Ndays = 350; % JBR - total number of days to calculate correction for (day-based approach)

