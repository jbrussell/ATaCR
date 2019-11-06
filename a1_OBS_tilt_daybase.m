% calculate the compliance and tilt effect for each station for each event.
% derive the tilt and compliance ndays (see parameter_tilt_evtbase.m) before the event.
% by PeiYing Patty Lin pylin.patty@gmail.com 201405
%
% JBR 2/22/18 : Modified for the purpose of removing tilt and compliance on
% the entire dataset for ambient noise applications. Instead of providing
% an "event" list, just calculate daily transfer functions.

clear;
% addpath ( '../EVTDATA_COMPLIANCE/');
parameter_tilt_evtbase;

dbpath = WORKINGdir; %'/Users/helenj/Cascadia/EARTHQUAKES/YEAR3/COMPLIANCESED/J42C/';
outpath = compliancematpath; %['/Users/jrussel/RESEARCH/PROJ_ENAM/TILTCOMP/OUT_DAYBASE/',PROJ,'/']; %'/Users/helenj/Cascadia/EARTHQUAKES/COMPLIANCE/';
figoutpath=['./',PROJ,'/figs/transfer/']; %'/Users/helenj/Cascadia/MATLABcodes/RF/Figures/COMPLIANCE/';


if ~exist(outpath)
    mkdir(outpath);
end
 
npts = T/dt;
Ppower = nextpow2(npts);
NFFT = 2^Ppower;
npad0 = (NFFT-npts)/2;
samprate = 1/dt; %common named Fs
f = samprate/2*linspace(0,1,NFFT/2+1);
Twin=1:(86400-T)/nwin:86400;
iptwin = (Twin-1)*samprate+1;

%datapath = [dbpath, char(eventids(ie)),'/']; % need to change
datapath = [dbpath];
%disp(datapath)
evtoutpath = [outpath];
%evtoutpath = [outpath,char(eventids(ie)),'/']; % need to change
if ~exist(evtoutpath)
    mkdir(evtoutpath);
end
if ~exist(figoutpath)
    mkdir(figoutpath);
end
for ista = 1: length(stalist) % 8:8  %12: 12 %length(Zsaclist)
    staname=char(stalist(ista,:));
    display(['Working on station : ',staname]);
    %     Zsaclist = dir([datapath,y,s,'*',chz]); 
%     Zsaclist = dir(sprintf('%s%s/%s/*%s.sac',datapath,staname,day,chz));
%     Z= readsac(fullfile(datapath,Zsaclist(ista).name)); 
%     staname = Z.KSTNM;
    
    list1 = dir([dbpath,staname,'/*',chz,'.sac']);
    for iday = 1 : length(list1)
        file1 = list1(iday).name;
        s = file1(6:22);        
%         s = num2str(iday,'%03d');
        if exist([evtoutpath,staname,'/',staname,'_',s,'_compliance.mat'])
            display(['skipping ',staname,' ',s]);
            continue
        else
            disp([staname,' ',s]);
        end


        % ---------- open zeros array for spectrum and conj(spectrum) for 4 channels -----
        spectrum_Z=zeros(nday*nwin,length(f));
        spectrum_H1=zeros(nday*nwin,length(f));
        spectrum_H2=zeros(nday*nwin,length(f));
        spectrum_P=zeros(nday*nwin,length(f));
        cspectrum_Z=zeros(nday*nwin,length(f));
        cspectrum_H1=zeros(nday*nwin,length(f));
        cspectrum_H2=zeros(nday*nwin,length(f));
        cspectrum_P=zeros(nday*nwin,length(f));
        
            
            % ------ Finding Start Times -----%
            
            
%             Zsaclist = dir(sprintf('%s%s/%s/*%s.sac',datapath,staname,day,chz));
%             filename = dir([datapath,y,s,'*.',staname,'.*',chz]);
%             filename = dir([datapath,staname,'/',s,'/*',chz,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',chz,'.sac']);
            if isempty(filename)
                display(['No data for ',chz,' ',staname,' day ',s]);
                continue
            end
            Z=readsac(fullfile(datapath,staname,filename.name));
            zhr=Z.NZHOUR;
            zmin=Z.NZMIN;
            zsec=Z.NZSEC;
            zmsec=Z.NZMSEC;
            zstart=zhr*3600+zmin*60+zsec+zmsec/1000; %start time in seconds from 00:00:00.00
                
%             filename = dir([datapath,y,s,'*.',staname,'.*',ch1]);
%             filename = dir([datapath,staname,'/',s,'/*',ch1,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',ch1,'.sac']);
            if isempty(filename)
                display(['No data for ',ch1,' ',staname,' day ',s]);
                continue
            end
            H1= readsac(fullfile(datapath,staname,filename.name));
            h1hr=H1.NZHOUR;
            h1min=H1.NZMIN;
            h1sec=H1.NZSEC;
            h1msec=H1.NZMSEC;
            h1start=h1hr*3600+h1min*60+h1sec+h1msec/1000;
            
%             filename = dir([datapath,y,s,'*.',staname,'.*',ch2]);
%             filename = dir([datapath,staname,'/',s,'/*',ch2,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',ch2,'.sac']);
            if isempty(filename)
                display(['No data for ',ch2,' ',staname,' day ',s]);
                continue
            end
            H2= readsac(fullfile(datapath,staname,filename.name));
            h2hr=H2.NZHOUR;
            h2min=H2.NZMIN;
            h2sec=H2.NZSEC;
            h2msec=H2.NZMSEC;
            h2start=h2hr*3600+h2min*60+h2sec+h2msec/1000;
            
%             filename = dir([datapath,y,s,'*.',staname,'.*',chp]);
%             filename = dir([datapath,staname,'/',s,'/*',chp,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',chp,'.sac']);
            if isempty(filename)
                display(['No data for ',chp,' ',staname,' day ',s]);
                continue
            end
            P= readsac(fullfile(datapath,staname,filename.name));
            phr=P.NZHOUR;
            pmin=P.NZMIN;
            psec=P.NZSEC;
            pmsec=P.NZMSEC;
            pstart=phr*3600+pmin*60+psec+pmsec/1000;
            
            
            %nptsZ = Z.NPTS;
            %nptsH1 = H1.NPTS;
            %nptsH2 = H2.NPTS;
            %nptsP = P.NPTS;
            
            % ---------- read DATA -------------
            
%             filename = dir([datapath,y,s,'*.',staname,'.*',chz]);
%             filename = dir([datapath,staname,'/',s,'/*',chz,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',chz,'.sac']);
           	[vec_tz,Zraw] = readsac(fullfile(datapath,staname,filename.name));
            
%             filename = dir([datapath,y,s,'*.',staname,'.*',ch1]);
%             filename = dir([datapath,staname,'/',s,'/*',ch1,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',ch1,'.sac']);
            [vec_t1,H1raw] = readsac(fullfile(datapath,staname,filename.name));
            
%             filename = dir([datapath,y,s,'*.',staname,'.*',ch2]);
%             filename = dir([datapath,staname,'/',s,'/*',ch2,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',ch2,'.sac']);
            [vec_t2,H2raw] = readsac(fullfile(datapath,staname,filename.name));
            
%             filename = dir([datapath,y,s,'*.',staname,'.*',chp]);
%             filename = dir([datapath,staname,'/',s,'/*',chp,'.sac']);
            filename = dir([datapath,staname,'/',staname,'.',s,'.','*',chp,'.sac']);
            [vec_tp,Praw] = readsac(fullfile(datapath,staname,filename.name));
            
            % ---------- aligning by start times -------------
            
            vec_tz=vec_tz+zstart;
            vec_t1=vec_t1+h1start;
            vec_t2=vec_t2+h2start;
            vec_tp=vec_tp+pstart;
            
            starttimes=[zstart,h1start,h2start,pstart]
            %npoints=[nptsZ,nptsH1,nptsH2,nptsP];
            [sorttimes,I]=sort(starttimes);
            
            if I(4)==1
                inttime=vec_tz;
            elseif I(4)==2
                inttime=vec_t1;
            elseif I(4)==3
                inttime=vec_t2;
            elseif I(4)==4
                inttime=vec_tp;
            end

%             Zraw  =  detrend(Zraw,'constant');
%             Zraw  =  detrend(Zraw,'linear');
%             H1raw =  detrend(H1raw,'constant');
%             H1raw =  detrend(H1raw,'linear');
%             H2raw =  detrend(H2raw,'constant');
%             H2raw =  detrend(H2raw,'linear');
%             Praw  =  detrend(Praw,'constant');
%             Praw  =  detrend(Praw,'linear');
            
            Zrawint = interp1(vec_tz,Zraw,inttime);
            H1rawint = interp1(vec_t1,H1raw,inttime);
            H2rawint = interp1(vec_t2,H2raw,inttime);
            Prawint = interp1(vec_tp,Praw,inttime);
            
            Zrawint(isnan(Zrawint)) = 0;
            H1rawint(isnan(H1rawint)) = 0;
            H2rawint(isnan(H2rawint)) = 0;
            Prawint(isnan(Prawint)) = 0;
            
            %Plotting to check that data has not been redacted and times
            %are aligned correctly
%             figure(1)
%             clf
%             subplot(411)
%             plot(vec_tz,cos_taper(Zraw),'-k');
%             hold on
%             plot(inttime,cos_taper(Zrawint),'-r');
%             xlabel('Time s');
%             Title(['Z Component,' char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname]);
%             subplot(412)
%             plot(vec_t1,cos_taper(H1raw),'-k');
%             hold on
%             plot(inttime,cos_taper(H1rawint),'-r');
%             xlabel('Time s');
%             Title(['H1 Component,' char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname]);
%             subplot(413)
%             plot(vec_t2,cos_taper(H2raw),'-k');
%             hold on
%             plot(inttime,cos_taper(H2rawint),'-r');
%             xlabel('Time s');
%             Title(['H2 Component,' char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname]);
%             subplot(414)
%             plot(vec_tp,cos_taper(Praw),'-k');
%             hold on
%             plot(inttime,cos_taper(Prawint),'-r');
%             xlabel('Time s');
%             Title(['P Component,' char(eventids(ie)),'_',char(ioridEVT(ie)),'-',staname]);
            
            %pause
     
            
            % ---------- get single-sided spectrum and conj(spectrum) 
            %            in each time window for 4 channels ----------
		

            for iwin = 1:length(iptwin) %iwin = (iiday-1)*nwin+1:iiday*nwin %Here might need to deal with problem of redacted data
%                 iiwin = iwin-(iiday-1)*nwin;
%                 pts_begin = iptwin(iiwin);
                pts_begin = iptwin(iwin);
                pts_end = pts_begin+npts-1;
                % Check window length
                if pts_begin > length(Zrawint) || pts_end > length(Zrawint)           
					disp('(Z) Points greater than the data... fixing window');
					pts_begin = length(Zrawint)-npts;
                    pts_end = pts_begin+npts;
                    %continue
                end 

                amp_Z  = Zrawint(pts_begin:pts_end);
                amp_Z  = detrend(amp_Z,'constant');
                amp_Z  = detrend(amp_Z,'linear');
                amp_Z  = cos_taper(amp_Z);
                amp_Z  = padarray(amp_Z,[0 npad0],'both');
                spectrum = fft(amp_Z,NFFT); 
                spectrum_Z(iwin,1:length(f)) = spectrum(1:NFFT/2+1);
                cspectrum_Z(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1));
                %loglog(f,2*abs(spectrum(1:NFFT/2+1))) 
                %hold on;
                %title('Single-Sided Amplitude Spectrum of y(t)')

                amp_H1 = H1rawint(pts_begin:pts_end);
                amp_H1 = detrend(amp_H1,'constant');
                amp_H1 = detrend(amp_H1,'linear');
                amp_H1 = cos_taper(amp_H1);
                amp_H1  = padarray(amp_H1,[0 npad0],'both');
                spectrum = fft(amp_H1,NFFT);
                spectrum_H1(iwin,1:length(f)) = spectrum(1:NFFT/2+1);
                cspectrum_H1(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1));
    %			%% --- quick plot check abs(spectrum)^2 == spectrum*cspectrum 
    %			loglog(f,(abs(spectrum(1:NFFT/2+1))).^2) 
    %			pause
    %			hold on;
    %			loglog(f,spectrum(1:NFFT/2+1) .* conj(spectrum(1:NFFT/2+1)),'r-')
    %			pause
    %			hold on;
    %			%% --- end quick plot

                amp_H2 = H2rawint(pts_begin:pts_end);
                amp_H2 = detrend(amp_H2,'constant');
                amp_H2 = detrend(amp_H2,'linear');
                amp_H2 = cos_taper(amp_H2);
                amp_H2  = padarray(amp_H2,[0 npad0],'both');
                spectrum = fft(amp_H2,NFFT);
                spectrum_H2(iwin,1:length(f)) = spectrum(1:NFFT/2+1);
                cspectrum_H2(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1));
                %loglog(f,2*abs(spectrum(1:NFFT/2+1))) 
                %hold on;
        

                amp_P  = Prawint(pts_begin:pts_end);
                amp_P  = detrend(amp_P,'constant');
                amp_P  = detrend(amp_P,'linear');
                amp_P  = cos_taper(amp_P);
                amp_P  = padarray(amp_P,[0 npad0],'both');
                spectrum = fft(amp_P,NFFT);
                spectrum_P(iwin,1:length(f)) = spectrum(1:NFFT/2+1); 
                cspectrum_P(iwin,1:length(f)) = conj(spectrum(1:NFFT/2+1));
                %loglog(f,2*abs(spectrum(1:NFFT/2+1))) 
                %hold on;
                
                
            end % for iwin
        
        % ---------- get auto-spectra g11,g22,gPP, gZZ 
        %            and cross-spectra g1Z,	g2Z,gPZ, g12, g1P, g2P 
        %            similar to eq(6) in Crawford and Webb, BSSA2000 ----------
        g11 = spectrum_H1.*cspectrum_H1;
        g11 = 2/nwin/(NFFT*dt)*sum(g11,1);

        g22 = spectrum_H2.*cspectrum_H2;
        g22 = 2/nwin/(NFFT*dt)*sum(g22,1);

        gPP = spectrum_P.*cspectrum_P;
        gPP = 2/nwin/(NFFT*dt)*sum(gPP,1);

        gZZ = spectrum_Z.*cspectrum_Z;
        gZZ = 2/nwin/(NFFT*dt)*sum(gZZ,1);

        g1Z = cspectrum_H1 .* spectrum_Z;
        g1Z = 2/nwin/(NFFT*dt)*sum(g1Z,1);

        g2Z = cspectrum_H2 .* spectrum_Z;
        g2Z = 2/nwin/(NFFT*dt)*sum(g2Z,1);    

        gPZ = cspectrum_P .* spectrum_Z;
        gPZ = 2/nwin/(NFFT*dt)*sum(gPZ,1);    

        g12 = cspectrum_H1 .* spectrum_H2;
        g12 = 2/nwin/(NFFT*dt)*sum(g12,1);

        g1P = cspectrum_H1 .* spectrum_P;
        g1P = 2/nwin/(NFFT*dt)*sum(g1P,1);

        g2P = cspectrum_H2 .* spectrum_P;
        g2P = 2/nwin/(NFFT*dt)*sum(g2P,1);
        
                     %plot coherences before removal
            
            cohhj_z1=abs(g1Z).^2./(g11.*gZZ);
            cohhj_z2=abs(g2Z).^2./(g22.*gZZ);
            cohhj_zp=abs(gPZ).^2./(gPP.*gZZ);
            cohhj_12=abs(g12).^2./(g11.*g22);
            cohhj_1p=abs(g1P).^2./(g11.*gPP);
            cohhj_2p=abs(g2P).^2./(g22.*gPP);
            figure(88); clf;
            ax1 = subplot(2,1,1);
            semilogx(f,smooth(abs(cohhj_z1),100),'r'); hold on
            semilogx(f,smooth(abs(cohhj_z2),100),'b'); hold on
            semilogx(f,smooth(abs(cohhj_zp),100),'k')
            xlim([1/100,1/2])
            ylim([0 1]);
            title('Coherence before removal');            
            legend('Z1','Z2','ZP','location','southeast');
            
            ax2 = subplot(2,1,2);
            semilogx(f,smooth(abs(cohhj_12),100),'r'); hold on
            semilogx(f,smooth(abs(cohhj_1p),100),'b'); hold on
            semilogx(f,smooth(abs(cohhj_2p),100),'k')
            legend('12','1p','2p','location','southeast');
            xlim([1/100,1/2])
            ylim([0 1]);
            xlabel('Frequency (Hz)');
            linkaxes([ax1,ax2],'xy');
            
            if ~exist([figoutpath,staname,'/'])
                mkdir([figoutpath,staname,'/'])
            end
            compliance_PS = [figoutpath,staname,'/',staname,'_',s,'_coherence_smooth.pdf'];
% %         print('-dpsc2',compliance_PS);
            save2pdf([compliance_PS],88,100);
        
        
        % ----------get the ZZ auto-spectrum after removing channel 1 gZZ_1, 
        %           then 2 gZZ_12, and then P, gZZ_12P
        %     plus  coherence between P and Z after removing effects of 1 and 2, gamPZ_12
        %     plus  ransfer function, l1Z,l12,l1P,l2P_1,l2Z_1,lPZ_12 
        %                             l1Z == conj(AZ1) in eq(8) in Crawford and Webb, BSSA2000 ----------
%         [gZZ_12P,gZZ_12,gZZ_1,gamPZ_12,l1Z,l12,l1P,l2P_1,l2Z_1,lPZ_12]=multicoher(g11,g22,gPP,gZZ,g1Z,g2Z,gPZ,g12,g1P,g2P,f);
        [gZZ_12P,gZZ_12,gZZ_1,gamPZ_12,gamPZ_1,gam2Z_1,l1Z,l12,l1P,l2P_1,l2Z_1,lPZ_12]=multicoher2(g11,g22,gPP,gZZ,g1Z,g2Z,gPZ,g12,g1P,g2P,f);

        

        % ---------- plot PS files for auto-Spectrum for gZZ, gZZ_12, and gZZ_12P ---------
        figure(102)
        clf
        ax1 = subplot(2,1,1);
        loglog(f,abs(g11),'b-');hold on;
        loglog(f,abs(g22),'r-');hold on;
        loglog(f,abs(gZZ),'k-');hold on;
        legend('g11','g22','gzz','location','southeast')
        title([s,'-',staname,'-auto-Spectrum for ZZ nwin ',num2str(nwin),' in a day, in ', num2str(nday),' days']);
        xlim([1/100,1/2])
        
        ax2 = subplot(2,1,2);
        loglog(f,abs(gZZ),'k-');hold on;
        loglog(f,abs(gZZ_12),'r-');hold on;
        loglog(f,abs(gZZ_12P),'b-');hold on;
        legend('gzz','gZZ-12','gZZ-12P','location','southeast')        
        title(['Z spectrum Corrected for tilt and compliance']);
        xlim([1/100,1/2])
        xlabel('Frequency (Hz)');  
        linkaxes([ax1,ax2],'xy');
        AutoSpectrumZZ_PS = [figoutpath,staname,'/',staname,'_',s,'_AutoSpectrumZZ.pdf'];
%         print('-dpsc2',AutoSpectrumZZ_PS);
        save2pdf([AutoSpectrumZZ_PS],102,100);

        
        figure(105); clf;
        subplot(2,1,1)
		loglog(f,smooth(abs(lPZ_12),100),'b-');
        xlim([1/100,1/2])
        xlabel('Frequency (Hz)');
        title([s,'-',staname,'-compliance: transfer function between P and Z after removing the effects of channels 1 and 2 '] )
        
        subplot(2,1,2)
        semilogx(f,smooth(cohhj_zp,100),'r'); hold on;
        semilogx(f,smooth(gamPZ_1,100),'b');
        semilogx(f,smooth(gamPZ_12,100),'--k');
        legend('ZP','ZP-2','ZP-1-2','location','southeast');
        title('Coherence between P and Z');
        xlim([1/100,1/2])
        ylim([0,1]);
        xlabel('Frequency (Hz)');        
        compliance_PS = [figoutpath,staname,'/',staname,'_',s,'_compliance.pdf'];
%         print('-dpsc2',compliance_PS);
        save2pdf([compliance_PS],105,100);
        
        
        figure(107); clf;
        subplot(2,1,1)
		semilogx(f,smooth(cohhj_z1,100),'r'); hold on
        semilogx(f,smooth(cohhj_z2,100),'b'); hold on
        semilogx(f,smooth(cohhj_zp,100),'k');
        semilogx(f,smooth(gamPZ_12,100),'--k');
        legend('Z1','Z2','ZP','ZP-1-2','location','southeast');
        xlim([1/100,1/2])
        ylim([0,1])
        xlabel('Frequency (Hz)');
        title([s,'-',staname,'Coherence before corrections'] )
        compliance_PS = [figoutpath,staname,'/',staname,'_',s,'_coherence_corr.pdf'];
%         print('-dpsc2',compliance_PS);
        save2pdf([compliance_PS],107,100);
        
        
        % ---------------- SAVE MATFILE ----------------
        if ~exist([evtoutpath,staname,'/'])
            mkdir([evtoutpath,staname,'/'])
        end
        compliance_MAT = [evtoutpath,staname,'/',staname,'_',s,'_compliance.mat'];
        save(compliance_MAT,'staname','f','l1Z','l12','l1P','l2P_1','l2Z_1','lPZ_12','gamPZ_12');
    end % for iday
    
   

end % end of loop ista

