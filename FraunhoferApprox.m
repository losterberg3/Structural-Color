% the code below is the derivation of the model. It also generates figures 1b-d of the paper

% grating pitch model (not included in paper due to Raman-Nath regime limitations)

lambda=0.400:0.001:0.430; % wavelengths

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4)); % refractive index
gp=[0.8;0.9;1;1.1]; % grating pitch
w1=0.1; %grating width
w2=0.2;
D1=w1./gp; % grating period
D2=w2./gp;

%navg1=(((gp-w1)./gp)*1)+(n_lambda.*(w1./gp)); 
%navg2=(((gp-w2)./gp)*1)+(n_lambda.*(w2./gp));

d1=200/1000; %thickness of each layer
d2=700/1000;

% fraunhofer approximation

dphi1=2*pi./lambda.*d1.*(n_lambda-1); % phase change for each layer
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

m=0;
eta1_0gp=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1)); % zeroth order diffraction for layer 1 and 2
eta2_0gp=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));

m=1;
eta1_1gp=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2); % first order diffraction for layer 1 and 2
eta2_1gp=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

eta_totalgp=eta2_0gp; 

figure
hold on
for i=1:length(gp)
    plot(1000*lambda,eta_totalgp(i,:),'Color',[0.85-(0.2*i) 0.85-(0.2*i) 1],'LineWidth',2)
end
hold off
ylim([0 1])
box on
xlabel('Wavelength (nm)','FontSize',16)
ylabel('Intensity (AU)','FontSize',16)
title 'Grating Pitch Model'
legend('800 nm','900 nm','1000 nm','1100 nm')
%% 

% height model for thinner and thicker layer (Figures 1b-c of paper)

lambda=0.400:0.001:0.700;

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp=1.1;
w1=0.1;
w2=0.2;
D1=w1./gp;
D2=w2./gp;

d1=[150;200;250;300]/1000; % thickness of each layer
d2=[500;550;600;650;700]/1000;

dphi1=2*pi./lambda.*d1.*(n_lambda-1);
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

m=0;
eta1_0h=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1)); %zeroth order diffraction for layer 1 and 2
eta2_0h=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));

m=1;
eta1_1h=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2); %first order diffraction for layer 1 and 2
eta2_1h=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

eta_totalh=eta2_0h;

figure
for i=1:length(d1)
    hold on
    plot(1000*lambda,eta1_0h(i,:),'Color',[0.85-(0.2*i) 0.85-(0.2*i) 1],'LineWidth',2)
end
hold off
ylim([0 1])
xlabel('Wavelength (nm)','FontSize',16)
ylabel('Intensity (AU)','FontSize',16)
title 'Figure 1b'
box on
legend('150 nm','200 nm','250 nm','300 nm')

figure
for i=1:length(d2)
    hold on
    plot(1000*lambda,eta_totalh(i,:),'Color',[1.1-(0.15*i) 1.1-(0.15*i) 1],'LineWidth',2)
end
hold off
ylim([0 1])
xlabel('Wavelength (nm)','FontSize',16)
ylabel('Intensity (AU)','FontSize',16)
title 'Figure 1c'
box on
legend('500 nm','550 nm','600 nm','650 nm','700 nm')


%% 
% height model using experimental heights (not included in paper)

lambda=0.400:0.001:0.700; % wavelengths

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp=1.1;
w1=0.1;
w2=0.2;
D1=w1./gp;
D2=w2./gp;

d1=[186;283;175;321;207;178]/1000; % thickness of each layer
d2=[691;839;694;907;716;662]/1000;
d2=d2-d1;

dphi1=2*pi./lambda.*d1.*(n_lambda-1); % phase change for each layer using fraunhofer approximation
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

m=0;
eta1_0h=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1)); %zeroth order diffraction for layer 1 and 2
eta2_0h=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));

m=1;
eta1_1h=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2); %first order diffraction for layer 1 and 2
eta2_1h=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

eta_totalh=eta2_0h;

figure
for i=1:length(d2)
    hold on
    plot(1000*lambda,eta_totalh(i,:),'LineWidth',2)
end
hold off
ylim([0 1])
xlabel 'Wavelength (nm)'
ylabel 'Intensity (AU)'
title 'Height Model Using Experimental Height'
legend('690 nm','839 nm','694 nm','907 nm','716 nm','662 nm')
%% 

% layers model (Figure 1d of paper)

lambda=0.400:0.001:0.700; % wavelengths

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4)); % refractive index
gp=1.1;
w1=0.1; % grating widths
w2=0.2;
w3=0.2;
w4=0.2;
D1=w1./gp;
D2=w2./gp;
D3=w3./gp;
D4=w4./gp;

d1=200/1000; % thickness of each layer
d2=700/1000;
d3=700/1000;
d4=700/1000;

dphi1=2*pi./lambda.*d1.*(n_lambda-1); % phase change of each layer from Fraunhofer approximation
dphi2=2*pi./lambda.*d2.*(n_lambda-1);
dphi3=2*pi./lambda.*d3.*(n_lambda-1);
dphi4=2*pi./lambda.*d4.*(n_lambda-1);

m=0;
eta1_0l=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1)); %zeroth order diffraction for layer 1, 2, 3, and 4
eta2_0l=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));
eta3_0l=1-(2.*D3)+(2.*D3.^2)+(2.*D3.*(1-D3).*cos(dphi3));
eta4_0l=1-(2.*D4)+(2.*D4.^2)+(2.*D4.*(1-D4).*cos(dphi4));

m=1;
eta1_1l=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2); %first order diffraction for layer 1, 2, 3, and 4
eta2_1l=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);
eta3_1l=4/pi^2/m^2.*(sin(pi.*m.*D3).^2).*(sin(dphi3./2).^2);
eta4_1l=4/pi^2/m^2.*(sin(pi.*m.*D4).^2).*(sin(dphi4./2).^2);

eta_total1=eta1_0l;  
eta_total2=eta2_0l.*eta1_0l;
eta_total3=eta2_0l.*eta1_0l.*eta3_0l;
eta_total4=eta2_0l.*eta1_0l.*eta3_0l.*eta4_0l;

figure
hold on
plot(1000*lambda,eta_total1,'Color',[0.85 0.85 1],'LineWidth',2)
plot(1000*lambda,eta_total2,'Color',[0.65 0.65 1],'LineWidth',2)
plot(1000*lambda,eta_total3,'Color',[0.45 0.45 1],'LineWidth',2)
plot(1000*lambda,eta_total4,'Color',[0.25 0.25 1],'LineWidth',2)
hold off
ylim([0 1])
box on
xlabel('Wavelength (nm)','FontSize',16)
ylabel('Intensity (AU)','FontSize',16)
title('Figure 1d','FontSize',16)
legend('1 layer','2 layers','3 layers','4 layers')
