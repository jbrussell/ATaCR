% transfer.m
% usage:
% out = transfer(data, npts, samples_per_second, 'pz_file', low_corner, high_corner)
% user made functions called:
% read_sac_pole_zero
% generate_response
% bp_bu_co
%
% procedes by reading the pole-zero file, fft'ing the data, computing the frequency array,
% computing the frequency dependent response, deconvolving (frequency
% domain dividing) the response from the data, and returns the data to the
% time domain with an ifft.
%
% does not use the fourier and ifourier routines to save time by only
% computing the necessary fft and ifft.
% 
% uses the high corner and low corner in a butterworth filter to stabilize
% the transfer function


function out = transfer(data,npts,sps,pz_file,low_corner,high_corner)
    a=fft(data,npts);
    [zz,pp,constant] = read_sac_pole_zero(pz_file);
    ss=length(a);
%   freq=zeros(1,ss);
    b=zeros(1,ss);
%    tmp=zeros(1,ss);
%    out=zeros(1,ss);
    for j=1:ss
%        freq = sps/ss*(j-1);
        freq = (j) / (1/sps * npts);
        tmp=generate_response(zz,pp,constant,freq);
        b(j) = zdiv(a(j),tmp); 
    end
%    tmp=ifft(tmp);
%    bb=(real(tmp));
    out = ifft(b,'symmetric');
    out = bp_bu_co(out,low_corner,high_corner,sps,4,2);
    
    return
        
