function [signal,ppm,B0_error]=fit_CEST( folder_path )
%% T2
here=pwd;
cd(folder_path)
file_names=dir('*mM*.txt');
number_of_samples=length(file_names);
offset=importdata('ppm.txt');
%offset=offset(5:end); 

 for q=1:number_of_samples
    s= importdata(file_names(q).name);
    s=s./s(4); 
    %s=s(5:end);
    [~,i]=max(1-s);
    B0_error(q,1)=offset(i);
    ppm(q,:)=offset-offset(i);
    signal(q,:)=s;
    pH(q,1)=str2double(file_names(q).name(3:6));
%Fittedpars(q,:)=fit_MSME(TE,s./max(s))'; %#ok<*AGROW>
 end
 

%Results=array2table([pH,Fittedpars],'VariableNames',{'pH','Estimated_T1','Lower_Bound','Upper_Bound','Rsquared'});

 cd(here)
%% plot T2 Error bar
% figure()
% y=Fittedpars(:,1);
% ub=Fittedpars(:,3)-Fittedpars(:,1);
% lb=Fittedpars(:,1)-Fittedpars(:,2);
%  errorbar(pH,y,lb,ub,'-o','LineWidth',2,'MarkerSize',10);
%  xlabel('pH','FontSize',20');
%   ylabel('Estimated T2 (seconds)','FontSize',20');
 
  
end
 



 
 


