% remove the tilt and compliance effect for Z components in EVENT based.
% read in the compliance and tilt information  
% and apply on the event data for all stations. 

% by PeiYing Patty Lin pylin.patty@gmail.com 201405

% add the compatibility that 
%   (1)time window length of the event sac files does not have to same as T.  
%   (2)if INSTRUMENTdir = 'NONE', still can run without making PSD plots.  
% by PeiYing Patty Lin  201406
%
% JBR 2/23/18 - Modified to read in full day files and chop them into
% length npts and apply correction to each section.
%


clear all;
setup_parameters;

isfigure = 1; 
iscompliance = 1; % 1 for tilt & compliance; 0 for tilt only
% compliancematpath = ['/Users/jrussel/RESEARCH/PROJ_ENAM/TILTCOMP/OUT_DAYBASE/',PROJ,'/']; %'/Users/helenj/Cascadia/EARTHQUAKES/YEAR1/FIXED/COMPLIANCE/';
outpath = ['/data/irma6/jrussel/YoungPacificORCA/SAC_1Hz_Zcorr_tiltcomp/']; % corrected seismograms; %'/Users/helenj/Cascadia/EARTHQUAKES/YEAR1/FIXED/SAC_out_tipr/';
figoutpath= ['./',PROJ,'/figs/SAC_1Hz_Zcorr_tiltcomp/']; %'/Users/helenj/Cascadia/MATLABcodes/RF/Figures/';


if ~exist(outpath)
    mkdir(outpath);
end
if ~exist(figoutpath)
    mkdir(figoutpath);
end

evtoutpath = [outpath];
 
npts = T/dt;
Ppower = nextpow2(npts);
NFFT = 2^Ppower;
samprate = 1/dt; %common named Fs
f = samprate/2*linspace(0,1,NFFT/2+1);


% matfiles = dir([compliancematpath,'*_compliance.mat']);

for ista = 1: length(stalist) % 8:8  %12: 12 %length(Zsaclist)
    staname=char(stalist(ista,:));
    display(['Working on station : ',staname]);

    list1 = dir([WORKINGdir,staname,'/*',chz,'.sac']);
    for iday = 1 : length(list1) %Ndays
        file1 = list1(iday).name;
        s = file1(6:22);
        %s = num2str(iday,'%03d');
        %Zsaclist = dir([WORKINGdir,staname,'/',s,'/','*',chz,'*']);
        %Z= readsac(fullfile(WORKINGdir,staname,s,Zsaclist.name));
        Zsaclist = dir([WORKINGdir,staname,'/',staname,'.',s,'.','*',chz,'.sac']);
        Z=readsac(fullfile(WORKINGdir,staname,Zsaclist.name));
        if isempty(Z)
            display(['No data for ',staname,' ',s]);
            continue
        end
        if Z.KSTNM ~= staname
            error('Wrong station!');
        end
        
        test = dir([evtoutpath,staname,'/',staname,'.*.',s,'*.sac']);
        if size(test,1) > 0
            display(['skipping ',staname,' ',s]);
            continue
        elseif size(test,1) == 0
            disp([s,' ',staname]);
        end

        matfile = dir([compliancematpath,staname,'/',staname,'_',s,'_compliance.mat']);
        
        matfilecell=struct2cell(matfile);
        [c d]=size(matfilecell);
        if d<=0
            display(['No compliance found for ',staname,' ',s,' ... skipping']);
            continue
        end
        
        filename = Z.FILENAME;
        filename = strrep(Z.FILENAME,chz,ch1);
        H1= readsac(fullfile(WORKINGdir,staname,filename));
        filename = strrep(Z.FILENAME,chz,ch2);
        H2= readsac(fullfile(WORKINGdir,staname,filename));
        filename = strrep(Z.FILENAME,chz,chp);
        P=  readsac(fullfile(WORKINGdir,staname,filename));
            
        evtsac_nptsZ = Z.NPTS;
        evtsac_nptsH1 = H1.NPTS;
        evtsac_nptsH2 = H2.NPTS;
        evtsac_nptsP = P.NPTS;
        npts_fils=[evtsac_nptsZ, evtsac_nptsH1, evtsac_nptsH2, evtsac_nptsP];
        if isequal(evtsac_nptsZ, evtsac_nptsH1, evtsac_nptsH2, evtsac_nptsP)
            evtsac_npts = Z.NPTS;
        else
            [sortpts,I]=sort(npts);
            evtsac_npts = sortpts(1); %This assumes they all have the same start time, but differing end times and will just cut off the longer traces.
        end
        
        clear l1Z l12 l1P l2P_1 l2Z_1 lPZ_12
        load(fullfile(compliancematpath,staname,matfile.name)); %Might not be correct path name
        
        filename = Z.FILENAME;
        [vec_tz,Zraw]  = readsac(fullfile(WORKINGdir,staname,filename));
        
        filename = strrep(Z.FILENAME,chz,ch1);
        [vec_t1,H1raw] = readsac(fullfile(WORKINGdir,staname,filename));
        
        filename = strrep(Z.FILENAME,chz,ch2);
        [vec_t2,H2raw] = readsac(fullfile(WORKINGdir,staname,filename));
        
        filename = strrep(Z.FILENAME,chz,chp);
        [vec_tp,Praw]  = readsac(fullfile(WORKINGdir,staname,filename));
        
        
        
        pts_begin = 1;
       	pts_end = evtsac_npts;
        
%         Ppower_n = nextpow2(npts(1));
%         NFFT = 2^Ppower_n;
% 		npad0 = floor((NFFT-evtsac_npts)/2);

%         npad0 = floor((NFFT-evtsac_npts)/2);
%         npad0 = floor((NFFT-6001)/2);
        
        ZCorrectwin = zeros(evtsac_npts,1);
        ampZ_12_nopad0 = zeros(evtsac_npts,1);
        ampZ_12P_nopad0 = zeros(evtsac_npts,1);
        amp_Z_nopad0 = zeros(evtsac_npts,1);
        winiter = ceil(evtsac_npts/npts);
        for iwin = 1:winiter % Loop over windows windows of correct length
            NFFT = 2^Ppower;
            pts_begin = 1 + (iwin-1)*npts;
            pts_end = pts_begin + npts-1;
            if pts_end > evtsac_npts
                pts_end = evtsac_npts;
            end
            npts_win = pts_end-pts_begin+1;
            
            if npts_win==npts
                npad0 = floor((NFFT-npts)/2);
            else
                npad0 = floor((NFFT-npts_win)/2);
            end
            
       
            % ---------get spectrum for 4 channels ----------
            amp_Z_win  = Zraw(pts_begin:pts_end);
            amp_Z_win  = detrend(amp_Z_win,'constant');
            amp_Z_win  = detrend(amp_Z_win,'linear');
            amp_Z_win  = cos_taper(amp_Z_win);
            qeqeqe = amp_Z_win;	
            amp_Z_win  = padarray(amp_Z_win,[0 npad0],'both');
            clear spectrum
            spectrum = fft(amp_Z_win,NFFT);
            spectrumZ = spectrum(1:NFFT/2+1);
    %     %			%% --- quick check plots for Zraw and ifft(spectrum_Zraw)
    %                 figure(222)
    %                 clf    			
    % 				ZZ = real(ifft(2*spectrumZ,NFFT));
    %         		Z_nopad0  = ZZ(npad0+1:npad0+1+evtsac_npts);
    % 				plot(Z_nopad0,'r');hold on;
    %     			pause;
    %     			plot(qeqeqe(pts_begin:pts_end),'b');
    %     %	        %% ---- end quick check plots

            amp_H1 = H1raw(pts_begin:pts_end);
            amp_H1 = detrend(amp_H1,'constant');
            amp_H1 = detrend(amp_H1,'linear');
            amp_H1 = cos_taper(amp_H1);
            amp_H1  = padarray(amp_H1,[0 npad0],'both');
            clear spectrum
            spectrum = fft(amp_H1,NFFT);
            spectrumH1 = spectrum(1:NFFT/2+1);

            amp_H2 = H2raw(pts_begin:pts_end);
            amp_H2 = detrend(amp_H2,'constant');
            amp_H2 = detrend(amp_H2,'linear');
            amp_H2 = cos_taper(amp_H2);
            amp_H2  = padarray(amp_H2,[0 npad0],'both');
            clear spectrum
            spectrum = fft(amp_H2,NFFT);
            spectrumH2 = spectrum(1:NFFT/2+1);


            amp_P  = Praw(pts_begin:pts_end);
            amp_P  = detrend(amp_P,'constant');
            amp_P  = detrend(amp_P,'linear');
            amp_P  = cos_taper(amp_P);
            amp_P  = padarray(amp_P,[0 npad0],'both');
            clear spectrum
            spectrum = fft(amp_P,NFFT);
            spectrumP = spectrum(1:NFFT/2+1);




            % ---------- removing tilt effect, 1st: removing chanel 1 effect ---------
            spectrumZ_1 = spectrumZ - ((l1Z) .* spectrumH1);

            % ---------- removing tilt effect, 2nd: removing chanel 2 effect ---------
            spectrumZ_12 = spectrumZ_1 - l2Z_1 .* (spectrumH2 - l12 .* spectrumH1);

            % ---------- removing P effect ---------
            spectrumZ_12P = spectrumZ_12 - lPZ_12 .* ( spectrumP - l1P.*spectrumH1 ...
            - l2P_1 .* (spectrumH2 - l12 .* spectrumH1)   );    





        %			%% --- quick plot for spectrumZ, spectrumZ_12 and spectrumZ_12P
        %	        loglog(f,smooth(2*abs(spectrumZ(1:NFFT/2+1)),100),'k-'); hold on;
        %	        loglog(f,smooth(2*abs(spectrumZ_12(1:NFFT/2+1)),100),'r-');hold on;
        %	        loglog(f,smooth(2*abs(spectrumZ_12P(1:NFFT/2+1)),100),'b-') 
        %	        %% --- end quick plot


            % ---------- ifft back spectrumZ_12 and spectrumZ_12P to time domain ---------
            % ---------- get rid off the pad0 cells ----------
            if npts_win==npts
                ampZ_12_win = real(ifft(2*spectrumZ_12,NFFT));
                ampZ_12P_win = real(ifft(2*spectrumZ_12P,NFFT));
                amp_Z_nopad0_win  = amp_Z_win(npad0+1:npad0+npts);
                ampZ_12_nopad0_win  = ampZ_12_win(npad0+1:npad0+npts);
                ampZ_12P_nopad0_win  = ampZ_12P_win(npad0+1:npad0+npts);
            elseif npts_win<npts
                ampZ_12_win = real(ifft(2*spectrumZ_12,NFFT));
                ampZ_12P_win = real(ifft(2*spectrumZ_12P,NFFT)); 
                amp_Z_nopad0_win  = amp_Z_win(npad0+1:npad0+npts_win);
                ampZ_12_nopad0_win  = ampZ_12_win(npad0+1:npad0+npts_win);
                ampZ_12P_nopad0_win  = ampZ_12P_win(npad0+1:npad0+npts_win);
            end

            amp_Z_nopad0_win = detrend(amp_Z_nopad0_win,'constant');
            amp_Z_nopad0_win = detrend(amp_Z_nopad0_win,'linear');	
            amp_Z_nopad0_win = cos_taper(amp_Z_nopad0_win);
        %			%% --- quick compare with original Z sacfile
        %            ZCorrect = Zsac;
        %            ZCorrect.DATA1 = amp_Z_nopad0;
        %            ZCorrect.FILENAME = strrep(Zsac.FILENAME,chz,'BCZ'); 
        %            status = writesac(ZCorrect);
        %			%% --- end quick check

            ampZ_12_nopad0_win = detrend(ampZ_12_nopad0_win,'constant');
            ampZ_12_nopad0_win = detrend(ampZ_12_nopad0_win,'linear');	
            ampZ_12_nopad0_win = cos_taper(ampZ_12_nopad0_win);

            ampZ_12P_nopad0_win = detrend(ampZ_12P_nopad0_win,'constant');
            ampZ_12P_nopad0_win = detrend(ampZ_12P_nopad0_win,'linear');	
            ampZ_12P_nopad0_win = cos_taper(ampZ_12P_nopad0_win);							

            if iscompliance
                ZCorrectwin(pts_begin:pts_end) = ampZ_12P_nopad0_win;
                ampZ_12P_nopad0(pts_begin:pts_end) = ampZ_12P_nopad0_win;
                ampZ_12_nopad0(pts_begin:pts_end) = ampZ_12_nopad0_win;
                amp_Z_nopad0(pts_begin:pts_end) = amp_Z_nopad0_win;
            else 
                ZCorrectwin(pts_begin:pts_end) = ampZ_12_nopad0_win;
                ampZ_12_nopad0(pts_begin:pts_end) = ampZ_12_nopad0_win;
                amp_Z_nopad0(pts_begin:pts_end) = amp_Z_nopad0_win;
            end
            
        end % iwin
        if length(ZCorrectwin)~=length(Z.DATA1)
            error('Number of points don''t match up!');
        end
        ZCorrect = Z;
        ZCorrect.DATA1 = ZCorrectwin;
		ZCorrect.FILENAME = strrep(Z.FILENAME,chz,'CZ'); %Will need to check string length for antelope
		status1 = writesac(ZCorrect);
        if ~exist([evtoutpath,staname,'/'])
            mkdir([evtoutpath,staname,'/']);
        end
        [status2, result] = system(['mv -f ', ZCorrect.FILENAME, ' ', evtoutpath,staname,'/' ]);
        

        amp_P_nopad0  = Praw(1:evtsac_npts);
        amp_P_nopad0  = detrend(amp_P_nopad0,'constant');
        amp_P_nopad0  = detrend(amp_P_nopad0,'linear');
        amp_P_nopad0  = cos_taper(amp_P_nopad0);			

        % ---------- plot seismograms -------------------------------------
		if isfigure
			fn = 1/2/dt;
			T1 = 20; T2= 50;	
			[b,a]=butter(2,[1/fn/T2,1/fn/T1]);
            
            amp_Z_filt  = filtfilt(b,a,amp_Z_nopad0);
            ampZ_12_filt  = filtfilt(b,a,ampZ_12_nopad0);
            ampZ_12P_filt  = filtfilt(b,a,ampZ_12P_nopad0);
            amp_P_filt  = filtfilt(b,a,amp_P_nopad0);
            if iscompliance
                figure(104)
                clf
                t = 0:dt:length(amp_Z_filt)-1;
                subplot(4,1,1)
                plot(t,amp_Z_filt(1:evtsac_npts),'k-');hold on;
                title(['amp-Z-filt (',num2str(T1),'-',num2str(T2),' s) ',staname,' ',s]);
                subplot(4,1,2)
                plot(t,ampZ_12_filt(1:evtsac_npts),'r-');hold on;
                title(['ampZ-12-filt ',staname,' ',s]);
                subplot(4,1,3)
                plot(t,ampZ_12P_filt(1:evtsac_npts),'b-');hold on;
                title(['ampZ-12P-filt ',staname,' ',s]);
                subplot(4,1,4)
                plot(t,amp_P_filt(1:evtsac_npts),'k-');hold on;
                title(['amp-P-filt ',staname,' ',s]);
                xlim([0 t(end)]);
            else
                figure(104)
                clf
                t = 0:dt:length(amp_Z_filt)-1;
                subplot(3,1,1)
                plot(t,amp_Z_filt(1:evtsac_npts),'k-');hold on;
                title(['amp-Z-filt (',num2str(T1),'-',num2str(T2),' s) ',staname,' ',s]);
                subplot(3,1,2)
                plot(t,ampZ_12_filt(1:evtsac_npts),'r-');hold on;
                title(['ampZ-12-filt ',staname,' ',s]);
                subplot(3,1,3)
                plot(t,amp_P_filt(1:evtsac_npts),'k-');hold on;
                title(['amp-P-filt ',staname,' ',s]);
                xlim([0 t(end)]);
            end                
            
            if ~exist([figoutpath,staname,'/'])
                mkdir([figoutpath,staname,'/'])
            end
            CorrectedSeismogram_PS = [figoutpath,staname,'/',staname,'_',s,'_CorrectedSeis.pdf'];
%             print('-dpsc2',CorrectedSeismogram_PS);
            save2pdf([CorrectedSeismogram_PS],104,250);
		end	
			
        
        
        
        % ---------- amp_Z remove instrument response 
        %            and transfer to ACC in oder to plot PSD ----------

        % ---------- get poles zeros values from SAC_pole_Zero file for Z comp ----------
		if isfigure && ~strcmp(INSTRUMENTdir,'NONE')	%Will need to edit this if I want to look at instrument response
			
			SACPoleZero = dir([INSTRUMENTdir,'SAC_PZs_ZA*',staname,'*',chz,'*']);
			SACPoleZero_Z = [INSTRUMENTdir,SACPoleZero.name];
			[zzeros,ppoles,ggain] =  read_sac_pole_zero(SACPoleZero_Z);
		
			figure(103)
			clf
			DISP= rm_SACPZ(amp_Z,zzeros,ppoles,ggain,samprate);
            VEL= time_differentiate(DISP, samprate);
			ACC= time_differentiate(VEL, samprate);
            
			xdft = fft(ACC,NFFT);
			xdft = xdft(1:NFFT/2+1);
			psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
			psdx(2:end-1) = 2*psdx(2:end-1);
			psdx = smooth(psdx,100);
			freq = 0:samprate/NFFT:samprate/2;
			semilogx(freq,10*log10(psdx),'k-'); hold on;


			% ---------- ampZ_12 remove instrument response 
			%            and transfer to ACC in oder to plot PSD ---------- 
			DISP= rm_SACPZ(ampZ_12,zzeros,ppoles,ggain,samprate);
			VEL= time_differentiate(DISP, samprate);
			ACC= time_differentiate(VEL, samprate);
            
			xdft = fft(ACC,NFFT);
			xdft = xdft(1:NFFT/2+1);
			psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
			psdx(2:end-1) = 2*psdx(2:end-1);
			psdx = smooth(psdx,100);
			freq = 0:samprate/NFFT:samprate/2;
			semilogx(freq,10*log10(psdx),'r-'); hold on;

			% ---------- ampZ_12P remove instrument response 
			%            and transfer to ACC in oder to plot PSD ----------					
			DISP= rm_SACPZ(ampZ_12P,zzeros,ppoles,ggain,samprate);
			VEL= time_differentiate(DISP, samprate);
			ACC= time_differentiate(VEL, samprate);

			xdft = fft(ACC,NFFT);
			xdft = xdft(1:NFFT/2+1);
			psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
			psdx(2:end-1) = 2*psdx(2:end-1);
			psdx = smooth(psdx,100);
			freq = 0:samprate/NFFT:samprate/2;
			semilogx(freq,10*log10(psdx),'b-');  hold on;

			legend('rawZ','Z-1-2','Z-1-2-P')
			% ---------- plot HLNM Peterson, 1993 ----------
			[LL,HH,FF] = noise_models(100); 
			semilogx(FF,HH,'color',[0.9 0.9 0.9]); hold on; 
			semilogx(FF,LL,'color',[0.9 0.9 0.9]); hold on;
			xx=[FF' fliplr(FF')];
			yy=[HH' fliplr(LL')];
			fill(xx,yy,[0.9 0.9 0.9],'FaceAlpha', 0.5,'EdgeColor','none');
			%fill(xx,yy,[0.9 0.9 0.9],'EdgeColor','none');
			xlim([0.01 20]);	  
			grid on;set(gca,'layer','top');
			xlabel('Hz');
			ylabel('Power[DB], 10 \ast log[m^{2}/s^{4}/Hz]');
			title(sprintf ('STA: %s',staname ));
			%ZZ_PS = [OUTPUTdir,'ZZ_OBScorrection_',sta,'_',num2str(year),'_',s,'.ps'];
			%print('-dpsc2',ZZ_PS);
            
            
            
            
            
            % plot DISP after removin instrument response
            figure(105)
            clf
            fn = 1/2/dt;
			T1 = 10; T2= 50;	
			[b,a]=butter(2,[1/fn/T2,1/fn/T1]);

            DISP_Z    = rm_SACPZ(amp_Z_nopad0 ,zzeros,ppoles,ggain,samprate);
            DISPZ_12  = rm_SACPZ(ampZ_12_nopad0 ,zzeros,ppoles,ggain,samprate);
            DISPZ_12P = rm_SACPZ(ampZ_12P_nopad0 ,zzeros,ppoles,ggain,samprate);
            
			amp_Z_filt  = filtfilt(b,a,DISP_Z);
			ampZ_12_filt  = filtfilt(b,a,DISPZ_12);
			ampZ_12P_filt  = filtfilt(b,a,DISPZ_12P);



			subplot(4,1,1)
			plot(amp_Z_filt(1:evtsac_npts),'k-');hold on;
			title('DISP-Z-filt');
			subplot(4,1,2)
			plot(ampZ_12_filt(1:evtsac_npts),'r-');hold on;
			title('DISPZ-12-filt');
			subplot(4,1,3)
			plot(ampZ_12P_filt(1:evtsac_npts),'b-');hold on;
			title('DISPZ-12P-filt');
            
            subplot(4,1,4)
			plot(amp_Z_filt(1:evtsac_npts),'k-');hold on;
			plot(ampZ_12_filt(1:evtsac_npts),'r-');hold on;
			plot(ampZ_12P_filt(1:evtsac_npts),'b-');hold on;
            title('plot 3 traces together');
            
            else
            figure(103)
			clf
            NFFT = evtsac_npts;
			DISP = amp_Z_nopad0;
            VEL= time_differentiate(DISP, samprate);
			ACC= time_differentiate(VEL, samprate);
            
			xdft = fft(ACC,NFFT);
			xdft = xdft(1:NFFT/2+1);
			psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
			psdx(2:end-1) = 2*psdx(2:end-1);
			psdx = smooth(psdx,100);
			freq = 0:samprate/NFFT:samprate/2;
			semilogx(freq,10*log10(psdx),'k-'); hold on;


			%            and transfer to ACC in oder to plot PSD ---------- 
			DISP= ampZ_12_nopad0;
			VEL= time_differentiate(DISP, samprate);
			ACC= time_differentiate(VEL, samprate);
            
			xdft = fft(ACC,NFFT);
			xdft = xdft(1:NFFT/2+1);
			psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
			psdx(2:end-1) = 2*psdx(2:end-1);
			psdx = smooth(psdx,100);
			freq = 0:samprate/NFFT:samprate/2;
			semilogx(freq,10*log10(psdx),'r-'); hold on;

			%            and transfer to ACC in oder to plot PSD ----------					
			DISP= ampZ_12P_nopad0;
			VEL= time_differentiate(DISP, samprate);
			ACC= time_differentiate(VEL, samprate);

			xdft = fft(ACC,NFFT);
			xdft = xdft(1:NFFT/2+1);
			psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
			psdx(2:end-1) = 2*psdx(2:end-1);
			psdx = smooth(psdx,100);
			freq = 0:samprate/NFFT:samprate/2;
            if iscompliance
                semilogx(freq,10*log10(psdx),'b-');  hold on;
                legend('rawZ','Z-1-2','Z-1-2-P','location','southeast')
            else
                legend('rawZ','Z-1-2','location','southeast')
            end
			% ---------- plot HLNM Peterson, 1993 ----------
			[LL,HH,FF] = noise_models(100); 
			semilogx(FF,HH,'color',[0.9 0.9 0.9]); hold on; 
			semilogx(FF,LL,'color',[0.9 0.9 0.9]); hold on;
			xx=[FF' fliplr(FF')];
			yy=[HH' fliplr(LL')];
			fill(xx,yy,[0.9 0.9 0.9],'FaceAlpha', 0.5,'EdgeColor','none');
			%fill(xx,yy,[0.9 0.9 0.9],'EdgeColor','none');
			xlim([1/100 1/2]);	  
			grid on;set(gca,'layer','top');
			xlabel('Hz');
			ylabel('Power[DB], 10 \ast log[m^{2}/s^{4}/Hz]');
			title(sprintf ('STA: %s %s',staname,s ));
            
            ZZ_PS = [figoutpath,staname,'/',staname,'_',s,'_ZZ_OBScorrection.pdf'];
%             print('-dpsc2',CorrectedSeismogram_PS);
            save2pdf([ZZ_PS],103,250);
		
        end
    end %end ista



end % end of loop ie

