close all, clc

g = 9.82;
lambda = linspace(0,0.15);
lambda = linspace(0,0.35,10000);
rho = 997;
sigma = 0.07275;
H = linspace(1,1,length(lambda));


%Dispersion relation
c = sqrt( ( g* lambda /(2* pi) + 2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
disp('The limiting cases')
%pure capillary wave
c_cap = sqrt( (  2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
%Pure gravity wave
c_grav = sqrt( ( g* lambda /(2* pi) ) .*tanh( 2*pi* H./lambda ) );
%Plots

figure(1);
% plot(lambda,c,lambda,c_cap,'-.',lambda,c_grav,'-.','LineWidth',1)
plot(lambda*100,c*100,lambda*100,c_cap*100,'-.',lambda*100,c_grav*100,'-.','LineWidth',1.5)



ylabel('$c$ [cm/s]','Interpreter','latex')
xlabel('$\lambda$ [cm]','Interpreter','latex')

legend('$c$','$c_{capillary}$','$c_{gravity}$','Interpreter','latex')
%     ylim([0 100])
%     xlim([0 6])
        ylim([0 100])
    xlim([0 20])
% disp('Shallow water approx, kH<<1, leads to tanh(kH)~kH')
% disp('Deep water approx, kH>>1, leads to tanh(kH)~1')
% %Deep water surface waves ( kH>>1 )
% c_cap_deep = sqrt( (2* pi* sigma)./ (rho*lambda)); %Pure capillary wave  (dispersive)
% c_grav_deep = sqrt( (g* lambda)./ (2*pi) ); %Pure gravity wave (dispersive)
% %Shallow water surface waves ( kH<<1 )
% c_cap_shallow = 2* pi./lambda .*sqrt( sigma.*H/ rho); %Pure capillary wave (dispersive)
% c_grav_shallow = sqrt(g* H); % Pure gravity wave (non dispersive)
lambda_min = 2*pi*sqrt(sigma/(rho*g))*10^2
lambda_cap = sqrt(sigma/(rho*g))*10^2
%2*pi*sqrt(lambda_)


set(gca,'fontsize',15) 
%figName = sprintf('dispersion_relation.png');
%exportgraphics(gca,figName,'Resolution',600)
%% accuracy of deep water approximation
clear, close, clc
g = 9.82;
rho = 997;
sigma = 0.07275;
lambda = 0.07;
%lambda=4e-3;
k = 2*pi/lambda;
c = sqrt( ( g/k + k*sigma/ rho ));
c_grav = sqrt( g/k);
c_cap = sqrt(k*sigma/ rho) ;
accuracy = (1-c_grav/c)*100; %accuracy in procent
disp('accuracy in %')
disp(accuracy)
disp('Surface tension only have an effect for wavelengths <7cm, this holds with an accuracy of  ~3%.')


 
%% figure of two sinusodial waves with slightly different wavenumbers 
clear, close, clc
omega = 6;% rad/sec
domega = 0.1; %rad/sec
k = 6;% rad/meter
dk = 0.3;% rad/meter.
a = 0.5;
t=linspace(0,1,1000);
x = linspace(-15,15,(length(t)));
eta = 2*a*cos(0.5*(dk.*x-domega.*t)).*cos(k.*x-omega.*t);
figure(2)
plot(x,eta,'LineWidth',1.5)
hold on
plot(x,2*a*cos(0.5*(dk.*x-domega.*t)),':r',x,-2*a*cos(0.5*(dk.*x-domega.*t)),':r','LineWidth',1.5)
xlabel('$x$','Interpreter','latex')
ylabel('$\eta$','Interpreter','latex')
set(gca,'fontsize',15) 
%figName = sprintf('wavepacket.png');
%exportgraphics(gca,figName,'Resolution',600)



%% Hyperbolic functions
clear, close, clc
x = linspace(0,2.5);
y = linspace(1,1);
y1 = cosh(x);
y2 = sinh(x);
y3 = tanh(x);

 v=1*0.97; %3procent avvikelse från 1
 v=0.94138; %detta ger 3%avvikelse från 1, pga roten ur!
 [~,idx] = (min(abs(y3 - v))); %hitta värdet i tanh(kH) som är närmast att avvika 3% från 1

figure(1)
plot(x,y1,x,y2,x,y3,x,y,':','LineWidth',1.5)
hold on
plot(x(idx),y3(idx),'or','linewidth',1.5)
ylabel('$y$','Interpreter','latex')
xlabel('$kH$','Interpreter','latex')
%legend('cosh($kH$)','sinh($kH$)','tanh($kH$)','Interpreter','latex')
legend('cosh($kH$)','sinh($kH$)','tanh($kH$)','$y$ = 1','tanh(1.75)','Interpreter','latex')
%For deep water approximation, already at kH=1.75, tanh(kH)=0.94138, since the phase speed in eq. (\ref{..}) is squared, this constitutes an deviation of 3% 
%For shallow water approximation, at kH=0.44, sqrt(tanh(kH)/0.44)=3%. i.e.
%for kH<0.44 the shallow water approximation holds with an 3%accuracy.
set(gca,'fontsize',15) 


% figName = sprintf('hyperbolic_functions.png');
% exportgraphics(gca,figName,'Resolution',600)

%%
clear, close, clc
g = 9.82;
lambda = 7.5e-2;
rho = 997;
sigma = 0.07275;
H = linspace(1,1,length(lambda));

c_grav_deep = sqrt( (g* lambda)./ (2*pi) );
c_g_min=c_grav_deep/2

t=4;
v=40;
s=v*t
lambda=linspace(0,3);
H_deep=0.28*lambda % [cm]
%For lambda_max<20cm, 5.6cm<H_deep, H_shallow<1.4cm
H_shallow = 0.07*lambda %[cm]


%% Perturbation pressure dependence on depth
clear, close, clc
H=1;
g = 9.82;
a=0.01;
lambda = 0.05;
k=2*pi/lambda;
rho = 997;
z = -1*linspace(0,0.1,1000);
p_prim = rho* g* a* cosh(k* (z+H))./cosh(k*H);

% figure(1)
% plot(p_prim/101.3,z/lambda) %Plot without x or t dependence, i have skipped the cos() term.
% xlabel("$P'/P_{atm}$ [-]",'Interpreter','latex')
% ylabel('z/$\lambda$ [-]','Interpreter','latex')
% grid on
% grid minor
% set(gca,'fontsize',12) 

%plot(-z/lambda,p_prim/(101.3e3),'LineWidth',1)
plot(-z/lambda,p_prim/p_prim(1),'LineWidth',1.5)
ylabel("$p'/p_{atm}$ [-]",'Interpreter','latex')
ylabel("$p'/p'(\eta)$ [-]",'Interpreter','latex')
xlabel('$\mid z\mid / \lambda$ [-]','Interpreter','latex')
grid on
grid minor
set(gca,'fontsize',15) 
% figName = sprintf('perturbation_pressure.png');
% exportgraphics(gca,figName,'Resolution',600)
%% water depth conditions for H with varying lambda
clear, close, clc
%Measurements in cm
lambda = linspace(0,15,1000);
lambda_cap = linspace(1.7,1.7,length(lambda));
lambda_gravdomin = linspace(7,7,length(lambda)); %Gravity dominating at wavelengths lambda>7cm
lambda_surftendomin = linspace(0.4,0.4,length(lambda)); %Surface tension dominating at wavelengths lambda<0.4cm
H_deep = 0.28*lambda;
H_shallow = 0.07*lambda;
H = linspace(0,H_deep(end),length(lambda));
figure; plot(lambda,H_deep,lambda,H_shallow,lambda_gravdomin,H,'-.',lambda_surftendomin,H,'-.','LineWidth',1.5)
legend('Deep, $H>0.28\lambda$','Shallow, $H<0.07\lambda$','Gravity dominates $\lambda>7$ cm','Surface tension dominates $\lambda<0.4$ cm','Interpreter','latex')
ylabel('$H$ [cm]','Interpreter','latex')
xlabel('$\lambda$ [cm]','Interpreter','latex')
set(gca,'fontsize',15) 
ylim([H(1) H(end)])

% figName = sprintf('inequality_plot.png');
% exportgraphics(gca,figName,'Resolution',600)


%% Camera angle of view
clear, close, clc
f = 50; %focal length 50mm
FOV = 430 ; %FOV = tank diameter in [mm]
d_H = 17.92; %sensor width horisonal (H) 17.92mm
d_V = 14.34; %sensor width vertical (V) mm
d_D = 22.9; %sensor width diagonal mm
d = [d_H,d_V,d_D];
disp('AOV in degrees: Horisontal, Vertical, Diagonal')
AOV = 2*atan(d/(2*f))*180/pi %Angle of view measured in degrees
disp('Camera placement above surface if FOV = tank diameter')
H = FOV ./ (2*tand(AOV/2))


%% Wave frequency for different lambda, to decide how fast FPS on the camera
clear, close, clc
g = 9.82;
lambda = linspace(0,0.15,500);
rho = 997;
sigma = 0.07275;
H = linspace(1,1,length(lambda));

%Dispersion relation
c = sqrt( ( g* lambda /(2* pi) + 2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );

f = c./lambda; %frequency [Hz]

plot(lambda*100,f)
ylabel('$f$ [Hz]','Interpreter','latex')
xlabel('$\lambda$ [cm]','Interpreter','latex')
set(gca,'fontsize',12) 
xlim([0 10])
lambda_grav_cap=[7,0.4];
c_grav_cap = [34,35]; %this is just approximat from the plot
hold on
f_grav_cap = c_grav_cap./lambda_grav_cap
scatter(lambda_grav_cap,f_grav_cap,'filled')
ylim([0 150])
legend('Frequencies between markers affected by both gravity and surf.tension')