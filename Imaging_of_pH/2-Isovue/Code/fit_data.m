function []=fit_data( folder_path,processing_function)
%% T1
here=pwd;
cd(folder_path)
file_names=dir('pH*');
number_of_samples=length(file_names);
TR=importdata('TR_ms.txt')./1000;

 for q=1:number_of_samples
    s= importdata(file_names(q).name);
    pH(q,1)=str2double(file_names(q).name(3:6));
[FittingPars(:,q)] =  fit_variable_TR(TR,s./max(s));

 end
 cd(here)
%% plot T1 Error bar
figure()
y=FittingPars(1,:)';
ub=FittingPars(3,:)-FittingPars(1,:); ub=ub';
lb=FittingPars(1,:)-FittingPars(2,:); lb=lb';
 errorbar(pH,y,lb,ub,'-o','LineWidth',2,'MarkerSize',10);
 xlabel('pH','FontSize',20');
  ylabel('Estimated T1 (seconds)','FontSize',20');


end
 



 
 


