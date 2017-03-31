clc; close all; clear all;
here=pwd;
addpath(here);
%% T1
cd ../data/PBS_only/T1
file_names=dir('pH*');
number_of_samples=length(file_names);
TR=importdata('TR_ms.txt')./1000;
VTR_signal=zeros(length(TR),number_of_samples);
 T1func=@(pars,xdata) pars(1) .* (1-exp(-xdata./pars(2)));
 
 
 for q=1:number_of_samples
    s= importdata(file_names(q).name);
    pH(q,1)=str2num(file_names(q).name(3:6));
    VTR_signal(:,q)=s./max(s);  
mdl=fitnlm(TR,VTR_signal(:,q),T1func,[1,3]);
ci=coefCI(mdl); ci(2,:);
Estimated_T1(q,1)=mdl.Coefficients.Estimate(2);
Estimated_T1(q,2:3)=ci(2,:);
E(:,q)=lsqcurvefit(T1func,[1,3],TR,VTR_signal(:,q),[1,1],[1.1,6]);
 end
 cd(here)
%% plot T1 Error bar
figure()
y=Estimated_T1(:,1);
ub=Estimated_T1(:,3)-Estimated_T1(:,1);
lb=Estimated_T1(:,1)-Estimated_T1(:,2);
 errorbar(pH,y,lb,ub,'-o','LineWidth',2,'MarkerSize',10);
 xlabel('pH','FontSize',20');
  ylabel('Estimated T1 (seconds)','FontSize',20');
 
%%% ##########

%% T2
cd ../data/PBS_only/T2
file_names=dir('pH*');
number_of_samples=length(file_names);
TE=importdata('TE_ms.txt')./1000;
%MSME_signal=zeros(length(TR),number_of_samples);
 T2func=@(pars,xdata) pars(1) .*  exp(-xdata./pars(2)) ;
 
 
 for q=1:number_of_samples
    s= importdata(file_names(q).name);
    %pH(q,1)=str2double(file_names(q).name(3:6));
mdl=fitnlm(TE,s./max(s),T2func,[1,.4]);
ci=coefCI(mdl); ci(2,:);
Estimated_T2(q,1)=mdl.Coefficients.Estimate(2);
Estimated_T2(q,2:3)=ci(2,:);
 end
 cd(here)
%% plot T1 Error bar
figure()
y=Estimated_T2(:,1);
ub=Estimated_T2(:,3)-Estimated_T2(:,1);
lb=Estimated_T2(:,1)-Estimated_T2(:,2);
 errorbar(pH,y,lb,ub,'-o','LineWidth',2,'MarkerSize',10);
 xlabel('pH','FontSize',20');
  ylabel('Estimated T2 (seconds)','FontSize',20');
 %%
 rmpath(here)
 
 



 
 


