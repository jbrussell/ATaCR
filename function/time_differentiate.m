% time_differentiate.m
% like the diff command to compute the derivative in the time domain
% and the freq_differentiate.m functions, but has the same effect in the time domain
% usage:
% out=time_differentiate(in,sps)
% where sps is the number of samples per second
% c version:
% //time derivative (as sac 2 point algorithm)
% void time_derivative(double *in, double *out, double delta, int npts) {
%	int i;
%	for (i=0; i<npts-1;i++) {
%		out[i] = (in[i+1] - in[i]) / delta;
%	}
%	//last value is not filled by sac. Here I fill it with the same as the previous value to have a minimal effect on future usage
%	out[npts-1] = out[npts-2];
%	
%	return;
% }
% //END

function out = time_differentiate(in,sps)
	out = zeros(1,length(in));
	for i=1:length(in) - 1
		out(i) = (in(i+1) - in(i)) * sps;
	end
	out(length(in)) = out(length(in) - 1);
	
	return
