function [Lsum,L]=lorentzian(pars,ppm)
% Goal ==>  Calculate a Lorenztian function for an N number of pools.
% Sintaxis:
%           L=lorentzian(pars,ppm)
%
%% Inputs 
%
% A         = vector of amplitudes for each function, length=N
% W         = vector of widths for each function,length=N
% x0        = vector of centers for each function,length=N
% ppm       = vector of range over which the Lorentzian will be calculated;
%
%% The input pars is equal to of the following styles
% pars= [A,W,x0], where  rows= number of Lorentizians and cols=parameters
% pars= [AW,x0], where  rows= number of Lorentizians and cols=parameters
%% Author
% Julio Cardenas-Rodriguez
% cardenaj@email.arizona.edu
% The University of Arizona Cancer Center
% Tucson, AZ
% Version 1.0
% September, 2013.

%% Steps
%scalingfactor=pars(end);
%pars=pars(1:end-1);
% 1) Preallocate output
L=zeros(length(ppm),(length(pars))/3);

% 4) Loop through each parameter
pools=(length(pars))/3;

for n=1:1:pools
P=pars(n:pools:end);
fwhm=P(2).^2/4;
L(:,n)=  (P(1).*fwhm) ./   (  fwhm + (ppm-P(3)).^2   )   ;
end
%
%
Lsum=sum(L,2) ; 
end





