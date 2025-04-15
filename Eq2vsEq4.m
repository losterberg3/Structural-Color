% The code below compares the integral from equation 2 of the paper with
% the analytical expression of equation 4 of the paper and generates Figure
% 1a of the paper

% Fraunhofer integral with scalar theory (equation 2)

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

% Fraunhofer analytical expression

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
