function [Hpsd1]=get_psd_nodb(obs1,samprate)

% get_psd.m  to return the counts power spectral densities for single
% channel

% demean the data
obs1=obs1-mean(obs1);

%subplot(211)
L=length(obs1);
Fs=samprate;
NFFT=2.^nextpow2(L);
obs1=double(obs1); % pylin.patty 08/20/2013 int32 array doesn't work!
Y=fft(obs1,NFFT)/L;
f=Fs/2*linspace(0,1,NFFT/2);
%loglog(f,2*abs(Y(1:NFFT/2))) 
%title('Single-Sided Amplitude Spectrum of y(t)')
%xlabel('Frequency (Hz)')
%ylabel('|Y(f)|')

% now plot the phase spectra
%subplot(212)
%semilogx(f, angle(Y(1:NFFT/2)).*180./pi)


% % Run without any filter, comment by pylin.patty 08/22/2013 =========
% % if(0)
% % % now convert to acceleration--velocity will have 0 0 0 in zeros
% % 
% % zeros = [ 0 0  -0.33 -5.2 -8.8];
% % poles = [0 -6.28 -6.28 -0.051 -0.06 -316 -0.43 -450 -350 (-150 -240i) (-150+240i)];
% % gain = 1.43e08;
% % 
% % topcoeff = poly(zeros);
% % botcoeff = poly(poles);
% % 
% % fax_rad = f.*2.*pi;
% % 
% % [H]=freqs(topcoeff, botcoeff, fax_rad);
% % figure(3)
% % subplot(212)
% % semilogx(f, (angle(H)).*180./pi, 'r')
% % subplot(211)
% % semilogx(f, 20.* log10(abs(H).*gain), 'r')
% % 
% % H2=[H,fliplr(H)];
% % 
% % H2=H2';
% % lambda = 1;
% % H2(isnan(H2))=complex(0);
% % corresp = conj(H2).*Y./(conj(H2).*H2 + lambda);
% % %corresp(1)=complex(0);
% % % convert from velocity to acceleration
% % 
% % Z=ifft(corresp, NFFT);
% % figure
% % subplot(211)
% % Zp=real(Z(1:L)).*L.*     3.397856604294215e+12;
% % plot(Zp)
% % 
% % end  %if(0)
% % ===============================================================


%figure
% now, power spectrum
h=spectrum.welch;
h.OverlapPercent = 50;
h.SegmentLength = 4096;
%fprintf(1, 'Length is %d, %s, " ", %s\n',length(obs1),obs_name, time0);
Hpsd1=psd(h, obs1, 'Fs', samprate);
%semilogx(Hpsd1.Frequencies, 10.*log10(Hpsd1.Data));






