function [CESTfit] =  cf_Lorentzian(cest_mri,ppm,Method,Control)
%% CestMRI
% Estimate the CEST effect of an N number of pools by non-linear LSQ.
% This function assumes that a Z-spectrum can be modeled as the sumation of
% Lorentzian fucntions in the frequency domain.
%
%%  SYNTAX
% [CESTfit] =  cf_Lorentzian(cest_mri,ppm,Method,Control)
% 
%%  INPUTS
%
% cest_mri: 3D Matrix.
%            where x=rows of image
%                  y=columns of images
%                  z=frequency offset of the transmitter
% ppm     : row vector of frequency offsets in ppm
% 
% Method: Structure with the following fields.
%            Method.Npools= number of lorentzians functions;
%            Method.range = row vector with first and data points to be
%                           analyzed
%            Method.x0    =  3*Npools by 3*Npools matrix that contains
%                          Method.x0(:,1)=Lower Limits for NLSQ fitting
%                          Method.x0(:,2)=Initial guess for NLSQ fitting
%                          Method.x0(:,3)=Upper Limits for NLSQ fitting
%            each column of is composed of:
%            rows 1     to N:   Lorentzian amplitudes for pools one to N
%            rows N+1   to N*2: Lorentzian width for pools one to N
%            rows 2*N+1 to N*3: Lorentzian offsets for pools one to N
%%  OUTPUTS
%
% CESTfit   Structure with teh follwinf elements:
%
%   CESTfit.pars:       P    arameters obtained after NLSQ fitting.
%          .residuals:   Residual of the fitting
%          .rsq:         Adjusted R-square
%          .cfitsum:     Predicted Lorentzian.
%          .cfitall:     Matrix of predicted invidual Lorentzians
%          .Zspectrum:   Experimental Zspectrum (full length)
%          .ppmadj:      Adjusted offsets to center water
%% References:  
%   1)  Sheth VR, Li Y, Chen LQ, Howison CM, Flask CA, Pagel MD.
%       Measuring in vivo tumor pHe with CEST-FISP MRI. 
%       Magn. Reson. Med., 2012, 67:760?768.  PMID 22028287.
%
%   2)  Liu G, Li Y, Sheth VR, Pagel MD.
%       Imaging in vivo extracellular pH with a Single PARACEST MRI Contrast Agent.
%       Molecular Imaging, 2012, 11(1):47-57.  PMID 21651182.
%
%
%% Author
%   Julio Cárdenas-Rodríguez
%       The University of Arizona
%       cardenaj@email.arizona.edu
%       ver 1.0 March, 2014.

%gaussFilter = fspecial('gaussian',[3 3],1);
%cest_mri= imfilter(cest_mri,gaussFilter,'same');

%% 1) Define variables from method structure:

x0=Method.x0;
point1=Method.range(1); point2=Method.range(2);
%% 4) Define options for NLSQ Method
% Fitting Method
options = optimoptions('lsqcurvefit');  
% Max. number of iterations
options.MaxIter=100E3; 
% Max. number of func. evaluations
options.MaxFunEvals=200;                
% Tolerance for NLSQ method
options.TolX=1e-4;             options.TolFun=1e-4;
% Turn iterative display off
options.Display='off';
options.Algorithm='trust-region-reflective';

%%  5)  Define parameters needed to allocate maps
ppm=ppm(point1:point2);   

%% 6) Perform Lorentzian Line Fitting only for Q voxels

    %% 6.1)Calculate Z spectrum and invert
    Signal= cest_mri(point1:point2);
    Zspectrum=Signal./Control;
    Zspectrum=1-Zspectrum;
    
    %% 6.2) Find Maximum and update ppm  to center Zspectrum
    [~,i]=max(Zspectrum);
    ppmadj=ppm-ppm(i);
    
    %% 6.3) Extrapolate and fit Zspectrum to a N-lorentzian model
ppmlong=linspace(min(ppmadj),max(ppmadj),1000); ppmlong=ppmlong';
Zlong = csaps(ppmadj,Zspectrum,1,ppmlong);
%%
[beta,~,R,~,cfoutput,~,Jacobian]=lsqcurvefit(@lorentzian,x0(:,2),ppmlong,Zlong,x0(:,1),x0(:,3),options);
        [Lsum,L]=lorentzian(beta,ppmadj);
            CESTfit.Lsum=Lsum;
            CESTfit.cfitall=zeros(length(Zspectrum),Method.Npools);
            CESTfit.cfitall=L;
            CESTfit.Zspectrum=Zspectrum;
            CESTfit.ppmadj=ppmadj;
            CESTfit.S=Signal;
            CESTfit.rsq=rsq(Zspectrum,Lsum'); 
            CESTfit.residuals=R;
            CESTfit.pars=beta;
            CESTfit.cfoutput=cfoutput; 
            
%             [Ypred,delta] = nlpredci(@lorentzian,ppmadj,CESTfit.pars,CESTfit.residuals,'Jacobian',Jacobian);
%            CESTpred=1-Ypred;
% CESTci=[1-(Ypred+delta),1-(Ypred-delta)];
%            
%           CESTfit.parsCI=nlparci(beta,R,'Jacobian',Jacobian); 
           
end
