%%Model using Kogelnik cwvt/Bragg diffraction/thick film assumption
% rotation angle model

lambda=0.350:0.001:0.750; %everything in micrometers

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp1=1.07;
gp2=1.07;
w1=0.1;
w2=0.18;

n_ave1=(n_lambda+1)/2; %(((gp1-w1)/gp1)*1)+(n_lambda.*(w1/gp1));
n_ave2=(n_lambda+1)/2; %(((gp2-w2)/gp2)*1)+(n_lambda.*(w2/gp2));

ng1=sqrt(((n_lambda.^2)*w1/gp1)+(1-(w1/gp1)));
ng2=sqrt(((n_lambda.^2)*w2/gp2)+(1-(w2/gp2)));

delta_n1=(n_lambda-1)/2; %w1*(n_lambda-1)/gp1; %accounting for optical thickness of the grating
delta_n2=(n_lambda-1)/2; %w2*(n_lambda-1)/gp2;

d2=(903.68-381.66)/1000; %thickness of each layer
d1=381.66/1000;

%theta1=asin(lambda./ng1/2/gp1); %using Bragg condition for the incident angle
%theta2=asin(lambda./ng2/2/gp2);

theta1=0; %incident angles on each layer
theta2=asin((lambda.*(n_ave1.^-1)/gp1)-sin(theta1)); %first order from the first layer

eta_p1=((sin((pi.*delta_n1.*d1)./(lambda.*cos(theta1)).*cos(2.*theta1))).^2);
eta_s1=((sin((pi.*delta_n1.*d1)./(lambda.*cos(theta1)))).^2);
eta_unpolarized1=0.5.*eta_p1+0.5.*eta_s1;
eta_p2=((sin((pi.*delta_n2.*d2)./lambda./cos(theta2).*cos(2.*theta2))).^2);
eta_s2=((sin((pi.*delta_n2.*d2)./lambda./cos(theta2))).^2);
eta_unpolarized2=0.5.*eta_p2+0.5.*eta_s2;

eta_0=eta_s1.*eta_p2;
eta_45=eta_unpolarized1.*eta_unpolarized2;
eta_90=eta_p1.*eta_s2;

figure
hold on
plot(lambda,eta_0,'LineWidth',2)
plot(lambda,eta_45,'LineWidth',2)
plot(lambda,eta_90,'LineWidth',2)
xlabel 'wavelength (micrometers)'
ylabel 'intensity (a.u.)'
title 'rotation angle'
legend('0','45','90')
%% 

% grating pitch model

lambda=0.350:0.001:0.750;

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp1=[0.77;0.87;0.97;1.07];
gp2=[0.77;0.87;0.97;1.07];
w1=0.1;
w2=0.18;

n_ave1=(n_lambda+1)/2; %(((gp1-w1)./gp1)*1)+(n_lambda.*(w1./gp1));
n_ave2=(n_lambda+1)/2; %(((gp2-w2)./gp2)*1)+(n_lambda.*(w2./gp2));

ng1=sqrt(((n_lambda.^2)*w1./gp1)+(1-(w1./gp1)));
ng2=sqrt(((n_lambda.^2)*w2./gp2)+(1-(w2./gp2)));

delta_n1=(n_lambda-1)/2; %w1*(n_lambda-1)/gp1; %accounting for optical thickness of the grating
delta_n2=(n_lambda-1)/2; %w2*(n_lambda-1)/gp2;

d2=(903.68-381.66)/1000; %thickness of each layer
d1=381.66/1000;

theta1=0; %incident angles on each layer
theta2=asin(lambda.*(n_ave1.^-1)./gp1); %first order from the first layer

eta_s1gp=((sin((pi.*delta_n1.*d1)./(lambda.*cos(theta1)))).^2);
eta_p2gp=((sin((pi.*delta_n2.*d2)./lambda./cos(theta2).*cos(2.*theta2))).^2);

eta_gp=eta_s1gp.*eta_p2gp;

figure
hold on
for i=1:4
    plot(lambda,eta_gp(i,:),'LineWidth',2)
end
hold off
xlabel 'wavelength (micrometers)'
ylabel 'intensity (a.u.)'
title 'grating pitch'
legend('770','870','970','1070')

% height model

lambda=0.350:0.001:0.750;

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp1=1.07;
gp2=1.07;
w1=0.1;
w2=0.18;

n_ave1=(n_lambda+1)/2; %(((gp1-w1)/gp1)*1)+(n_lambda.*(w1/gp1));
n_ave2=(n_lambda+1)/2; %(((gp2-w2)/gp2)*1)+(n_lambda.*(w2/gp2));

ng1=sqrt(((n_lambda.^2)*w1./gp1)+(1-(w1./gp1)));
ng2=sqrt(((n_lambda.^2)*w2./gp2)+(1-(w2./gp2)));

delta_n1=(n_lambda-1)/2; %w1*(n_lambda-1)/gp1; %accounting for optical thickness of the grating
delta_n2=(n_lambda-1)/2; %w2*(n_lambda-1)/gp2;

d1=[203.22;281.87;174.84;324.84;188.46;169.18]/1000; %thickness of each layer
d2=[495.06;560.25;523.98;578.84;524.77;486.95]/1000;

theta1=0; %incident angles on each layer
theta2=asin(lambda.*(n_ave1.^-1)./gp1); %first order from the first layer

eta_s1d=((sin((pi.*delta_n1.*d1)./(lambda.*cos(theta1)))).^2);
eta_p2d=((sin((pi.*delta_n2.*d2)./lambda./cos(theta2).*cos(2.*theta2))).^2);

eta_d=eta_s1d.*eta_p2d;

figure
for i=1:6
    hold on
    plot(lambda,eta_d(i,:),'LineWidth',2)
end
hold off
xlabel 'wavelength (micrometers)'
ylabel 'intensity (a.u.)'
title 'height'
legend('698','841','698','903','713','656')
        
%% duty cycle model for thin film assumption (Raman Nath grating)

% grating pitch model

lambda=0.400:0.001:0.700;

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp=[0.77;0.87;0.97;1.07];
w1=0.1;
w2=0.18;
D1=w1./gp;
D2=w2./gp;

d1=0.322; %thickness of each layer
d2=0.907-0.322;
%d2=[0.675;0.71;0.73;0.79];

dphi1=2*pi./lambda.*d1.*(n_lambda-1);
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

% different orders

m=0;
eta1_0=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1));
eta2_0=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));

m=1;
eta1_1=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_1=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=-1;
eta1_neg1=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_neg1=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=2;
eta1_2=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_2=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=-2;
eta1_neg2=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_neg2=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

eta_total=eta2_0; %+eta2_1+eta2_neg1; %+eta2_2+eta2_neg2; 

figure
hold on
for i=1:4
    plot(lambda,eta_total(i,:),'LineWidth',2)
end
hold off
%ylim([0 1])
xlabel 'wavelength (micrometers)'
ylabel 'intensity (a.u.)'
title 'grating pitch'
legend('770','870','970','1070')

%% 

% height model

lambda=0.400:0.001:0.700;

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp=1.07;
w1=0.1;
w2=0.18;
D1=w1./gp;
D2=w2./gp;

d1=[186;283;175;321;207;178]/1000; %thickness of each layer
d2=[691;839;694;907;716;662]/1000;
d2=d2-d1;
%d2=(400:75:700)'/1000;

dphi1=2*pi./lambda.*d1.*(n_lambda-1);
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

% different orders

m=0;
eta1_0=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1));
eta2_0=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));

m=1;
eta1_1=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_1=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=-1;
eta1_neg1=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_neg1=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=2;
eta1_2=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_2=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=-2;
eta1_neg2=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_neg2=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

eta_total=eta2_0+eta2_1+eta2_neg1; %+eta2_2+eta2_neg2;

figure
for i=1:6
    hold on
    plot(lambda*1000,eta_total(i,:),'LineWidth',2)
end
hold off
ylim([0 1])
xlabel 'Wavelength (nanometers)'
ylabel 'Intensity (a.u.)'
title 'Height Variation Analytical Model'
legend('690','839','694','907','716','662')
%legend('400 nm','475 nm','550 nm','625 nm','700 nm','900')

%% 
% light polarization model

lambda=0.400:0.001:0.700;

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp=1.07;
w1=0.1;
w2=0.18;
D1=w1./gp;
D2=w2./gp;
d1=0.322;
d2=0.907-0.322;

%np sin thetam + ni sin theta i = -m lambda / period, assume thin
%diffracting layer so that n=n'=1 for air
%theta1=asin(-m.*lambda/gp);

dphi1=2*pi./lambda.*d1.*(n_lambda-1);
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

%0 deg

% different orders

m=0;
theta1=asin(-m.*lambda/gp);
eta1_0=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1));
eta2_0=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));

m=1;
theta1=asin(-m.*lambda/gp);
eta1_1=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_1=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

m=-1;
theta1=asin(-m.*lambda/gp);
eta1_neg1=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_neg1=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);

eta_total=eta1_0.*eta2_0;

plot(lambda,eta_total)

%% Klein parameter

%for height variation
lambda=0.400:0.001:0.700; %everything in micrometers

A=1.52266;
B=0.000733;
C=-0.0000482;
n_lambda=A+(B*(lambda.^-2))+(C*(lambda.^-4));
gp=1.07;
w1=0.1;
w2=0.18;

d1=[186;283;175;322;208;178]/1000; %thickness of each layer
d2=[691;839;694;907;716;662]/1000;
d2=d2-d1;

n_ave1=(((gp-w1)/gp)*1)+(n_lambda.*(w1/gp));
n_ave2=(((gp-w2)/gp)*1)+(n_lambda.*(w2/gp));

Q1=2*pi*d1.*lambda/gp/gp./n_ave1;
Q2=2*pi*d2.*lambda/gp/gp./n_ave1;

figure
for i=1:6
    hold on
    plot(lambda,Q1(i,:),'LineWidth',2)
end
plot(lambda,zeros(length(lambda))+1)
hold off
ylim([0 3])
xlabel 'wavelength (micrometers)'
ylabel 'Q Parameter'
title 'height, layer 1'
legend('690','839','694','907','716','662')

figure
for i=1:6
    hold on
    plot(lambda,Q2(i,:),'LineWidth',2)
end
plot(lambda,zeros(length(lambda))+1)
hold off
ylim([0 3])
xlabel 'wavelength (micrometers)'
ylabel 'Q Parameter'
title 'height, layer 2'
legend('690','839','694','907','716','662')

%for gp variation

gp=[0.77;0.87;0.97;1.07];
w1=0.1;
w2=0.18;

d1=0.322; %thickness of each layer
d2=0.907-0.322;

n_ave1=(((gp-w1)./gp)*1)+(n_lambda.*(w1./gp));
n_ave2=(((gp-w2)./gp)*1)+(n_lambda.*(w2./gp));

Q1=2*pi*d1.*lambda./gp./gp./n_ave1;
Q2=2*pi*d2.*lambda./gp./gp./n_ave1;

figure
hold on
for i=1:4
    plot(lambda,Q1(i,:),'LineWidth',2)
end
plot(lambda,zeros(length(lambda))+1)
hold off
ylim([0 3])
xlabel 'wavelength (micrometers)'
ylabel 'Q Parameter'
title 'grating pitch, layer 1'
legend('770','870','970','1070')

figure
hold on
for i=1:4
    plot(lambda,Q2(i,:),'LineWidth',2)
end
plot(lambda,zeros(length(lambda))+1)
hold off
ylim([0 3])
xlabel 'wavelength (micrometers)'
ylabel 'Q Parameter'
title 'grating pitch, layer 2'
legend('770','870','970','1070')

%% Integral with scalar theory

gp=1.1;
w=0.1;
d1=0.2;
d2=0.7;
rect= @(X) 1+(0.52*rectangularPulse(-w/2,w/2,X-(gp/2)));
m=0;

eta1=zeros(1,16);
for lambda=20:35
    Em= @(X) exp(-2*pi*1i*((m*X/gp)+(d1*rect(X)/(lambda/50))))/gp; 
    q=integral(Em,0,gp);
    eta1(lambda-19)=q*conj(q);
end

gp=1.1;
w=0.2;
d1=0.2;
d2=0.7;
rect= @(X) 1+(0.52*rectangularPulse(-w/2,w/2,X-(gp/2)));
m=0;

eta2=zeros(1,16);
for lambda=20:35
    Em= @(X) exp(-2*pi*1i*((m*X/gp)+(d2*rect(X)/(lambda/50))))/gp; 
    q=integral(Em,0,gp);
    eta2(lambda-19)=q*conj(q);
end

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

d1=0.200; %thickness of each layer
d2=0.700;

dphi1=2*pi./lambda.*d1.*(n_lambda-1);
dphi2=2*pi./lambda.*d2.*(n_lambda-1);

eta1_0=1-(2.*D1)+(2.*D1.^2)+(2.*D1.*(1-D1).*cos(dphi1));
eta2_0=1-(2.*D2)+(2.*D2.^2)+(2.*D2.*(1-D2).*cos(dphi2));
m=3;
eta1_m=4/pi^2/m^2.*(sin(pi.*m.*D1).^2).*(sin(dphi1./2).^2);
eta2_m=4/pi^2/m^2.*(sin(pi.*m.*D2).^2).*(sin(dphi2./2).^2);
%% 

figure
hold on
plot(0.4:0.02:0.7,eta1,'LineWidth', 3)
xlabel 'Wavelength (Micrometers)'
ylabel 'Intensity (a.u.)'
title 'Diffraction Efficiency'
plot(0.4:0.02:0.7,eta2,'LineWidth', 3)
plot(lambda,eta1_0,'LineWidth', 3)
plot(lambda,eta2_0,'LineWidth', 3)
ylim([0 1])
legend('Thin Layer Integral Method','Thick Layer Integral Method','Thin Layer Analytical Method (Eq 4)','layer 2 analytical (cos func)')
hold off

