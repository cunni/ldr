%%%%%%%%%%%%%%%%%%%%
% John P Cunningham
% 2013
%
% calculates log(1 + exp(z)) and its derivative without overflow or
% underflow
% original code by Byron Yu
% modified by John P Cunningham 2007
% and then again in 2013
%%%%%%%%%%%%%%%%%%%%

function [ hz , dhz ] = log_one_plus_exp(z)

    scale = 1;

    z = scale*z;
    outUL = (z>10/scale);
    outLL = (z<-10/scale);
    inLim = ~outUL & ~outLL;
    dhz = nan(size(z));
    % return derivative
    dhz(outUL) = 1;
    dhz(outLL) = exp(z(outLL));
    dhz(inLim) = exp(z(inLim))./(1+exp(z(inLim)));
    
    hz = nan(size(z));
    % return logOnePlusExp
    hz(outUL) = z(outUL);
    hz(outLL) = exp(z(outLL));
    hz(inLim) = log(1 + exp(z(inLim)));
    
end
