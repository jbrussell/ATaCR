% zdiv.m
% function to do complex division
% should be already available in matlab, but seems not to work so well
% usage:
% out=zdiv(a+bi,c+di);

function c=zdiv(a,b)
    denom = real(b) * real(b) + imag(b) * imag(b)+1.2e-32;
    c =  (real(a)*real(b)+imag(a)*imag(b)) / denom + i * (imag(a)*real(b)-real(a)*imag(b)) / denom;
   
    return
 
