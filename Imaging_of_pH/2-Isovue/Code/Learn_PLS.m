clc; close all; clear all;
 load moore
y = moore(:,6);              % Response
X0 = moore(:,1:5);           % Original predictors
X1 = X0+10*randn(size(X0));  % Correlated predictors
X = [X0,X1];

%   Use plsregress to perform PLS regression with the same number 
%   of components as predictors, then plot the percentage variance 
%   explained in the response as a function of the number of components:

[XL,yl,XS,YS,beta,PCTVAR] = plsregress(X,y,10);
yfit = [ones(size(X,1),1) X]*beta;

subplot(121);
plot(1:10,cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y');
subplot(122);
scatter(y,yfit,100,'filled'); lsline