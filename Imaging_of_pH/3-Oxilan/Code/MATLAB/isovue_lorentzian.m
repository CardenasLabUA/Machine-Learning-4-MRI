%% Center data via Lorentzian fitting
% Prescribe number of pools for lorentzian fit
method.Npools=3;
center_range=[-5,10];

% Prescribe initial guesses and lower and upper bounds for amplitudes,
% widths, and offsets of lorentzian peaks from the 3 pools
method.x0=[0.4, 1, 1.5;... Water amplitude - 1
    0, 0.3, 1;... Oxilan pool 1 amplitude - 2
    0, 0.2, 1;... Oxilan hydroxyls amplitude - 3
    0, 0.3, 1;... isovue 5.6 - 4
    0, 4, 10;... Water width - 5
    0, 4, 15;... Oxilan pool 1 width - 6 
    0, 2, 15;... Oxilan hydroxyls width - 7
    0, 2, 15;... isovue 5.6 - 8
    -0.5, 0, 0.5;... Water offset - 9
    3.5, 4.2, 4.7;... Oxilan pool 1 offset - 10
    0.5, 2.5, 3.5;... Oxilan hydroxyls offset - 11
    4.7, 5.6, 6];... isovue 5.6 - 12

xpred=linspace(center_range(1),center_range(2),length(ppm));
for j=1:size(c_signal,1)
    
    control=find(round(ppm(j,:))==center_range(1),1);
    method.range=[1,length(c_signal(j,:))];
    fit=cf_Lorentzian(c_signal(j,:),ppm(j,:),method,c_signal(j,control));
    fit.pars(9:end)=fit.pars(9:end)-fit.pars(9);
    c_centered(j,:)=lorentzian(fit.pars,xpred)';
    c_rsq(1,j)=fit.rsq;
    pars(j,1)=fit.pars(2);
    pars(j,2)=fit.pars(4);
    pars(j,3)=fit.pars(6);
    pars(j,4)=fit.pars(8);
end

%% Plotting
figure('Name','Centered Data');

pHs = pH(1:9:end);
cool=reshape(c_centered',length(c_xpred),9,8);
for q=1:8
    subplot(3,3,q); plot(c_xpred,1-cool(:,:,q));
    title(['(Centered Data) pH = ',num2str(pHs(q))]); ylim([0 1]);
    xlabel('Offset (ppm)')
    ylabel('Mz/Mo')
end
subplot(3,3,9); plot(c_xpred,1-c_centered');title('All');