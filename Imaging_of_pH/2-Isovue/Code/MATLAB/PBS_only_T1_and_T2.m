clear all; close all; clc;
% Move to data directory
here=pwd;
%% T1
% names
cd ../data/PBS_only/T1/
file_names=dir('pH*');
TRInfo=dir('*TR*'); TR=importdata(TRInfo.name); TR=TR'./1000;

% number of files
num_files=length(file_names);
% TR

for q=1:num_files
S=importdata(file_names(q).name); 
end
cd(here)

%%