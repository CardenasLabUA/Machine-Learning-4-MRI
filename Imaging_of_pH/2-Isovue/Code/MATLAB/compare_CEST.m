clear all; close all;
% pH6.25		pH6.77		pH7.19		pH8.00
% pH6.02		pH6.60		pH7.03		pH7.40

[signal,ppm]=fit_CEST( '../data/pH6.02/CEST/' ); 
plot(ppm,signal,'or');% xlim([ -10 20]);
hold all;

[signal,ppm]=fit_CEST( '../data/pH8.00/CEST/' ); 
plot(ppm,signal,'xk'); %xlim([ -10 20]);