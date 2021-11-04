function Podosome_wave_dynamics

% pre-define all the parameters
phi = cos(pi/4);
x0 = 2800;                        % unit nm
Vr = 100;                        % unit nm/s
alpha = 1.5;                    % unit 1
beta = 0.5;                     % unit 1/nm
css = 20;                           % unit uM
gamma = 1*0.2;
ks = 200000;                        % unit pN/nm
Kf = 3.5;
kc = 40;
eta = 0.0;                     % unit pN*s/nm change Kf will change the xss
%tau = 0.0;                      % unit s
Vp0 = 70;                     % unit nm/s
Vd = 50;                       % unit nm/s
F0 = 1000;                       % unit pN   decrease this value increases xss
Fp0 = 20000;                     % unit pN
tfinal = 3000;                 % unit s

% define the podosome connectivity and positions
Na = 10;                           % number of array
%Npod = Na*(Na+1);                       % podosome number
d0 = 1500;                         % default podosome distance unit 1 um
[podoconnect, xpod, ypod, Npod] = PodoConnectivity_hexagon(Na, d0);

% Equilibirum point
rfs = 1;
rcs = ks/(kc+ks);
Vpss = Vp0+beta*css;
L1ss = (Vpss-Vd)*Fp0*rcs/ks/Vpss;
Fmss = Kf/(Kf-gamma)*(F0+(alpha-gamma/Kf/rfs)*ks*L1ss/phi);
xss = (ks*L1ss/phi/rfs-Fmss)/Kf+x0;
Lss = phi*xss;

% parameter groups
a = rcs*Fp0/Vpss/ks*(ks/Kf/rfs+1);
b = rcs*Fp0*phi/(Kf*Vpss*ks);
Tmem = eta/ks;
c = 1-gamma/Kf;
e = (alpha-gamma/Kf/rfs)*ks/phi;

tau =  40;
Real = (-a*c+b*e-tau)/2/(a*tau+beta*b);    % because denominator a-c<0
Imagine = sqrt(-4*c*(beta*b+tau*a)+(a*c-b*e+tau)^2)/2/(a*tau+beta*b);  
Period = 2*pi/abs(Imagine);  %~2*pi*sqrt(2*a*tau)

% main equations for six variables:  Fm, x
% Use an inconsistent initial condition to test initialization.
y0 = L1ss*ones(3*Npod,1);
y0(2:3:end)= Fmss;%-500;  
y0(3:3:end) = css; 
Midpod = 1; Midpodneigh = podoconnect{Midpod};
% for radial waves
% y0(Midpod*3-1)=Fmss-1500; y0(Midpodneigh*3-1)=Fmss;
% for random waves
y0(2:3:end)=Fmss-250*(rand(Npod,1)*2-1); 

BoundPodo = Npod-(Na-1)*6+1:Npod;
y0(3*BoundPodo-1)=Fmss; 



spara = [Kf, alpha, beta, ks, tau, Vp0, Vd, F0, Fp0, phi, gamma, rcs, rfs, d0, Npod];

% Use options to defineteh relative toleration
options = odeset('RelTol',1e-3);
tspan = 0:5:tfinal;
[t,y] = ode45(@(t,y) odefcn(t,y,spara, podoconnect), tspan, y0);
xcable = (ks*y(:,1:3:end)/phi/rfs-y(:,2:3:end))/Kf+x0;

% generate the movie
%movieplot2(Na, d0, xcable(end-300:1:end,:)*phi)

% plot the results
figure(1)
subplot(2,2,2)
plot(t, y(:,Midpod*3-2))%+xcable*phi)
hold on
ylabel('L1')
subplot(2,2,3)
plot(t, y(:,Midpod*3-1))
hold on
plot(t, ks*y(:,Midpod*3-2))
hold on
ylabel('Fm')
subplot(2,2,1)
plot(t, xcable*phi)
hold on
%plot(t, xss*ones(length(t),1))
ylabel('x')
subplot(2,2,4)
plot(t, y(:,Midpod*3))
hold on
%ylabel('Fp')
%plot(t, y(:,5)./y(:,3))

figure(2)
surf(xcable*phi)


end
%
% function definition
%
function dydt = odefcn(t,y,spara, podoconnect)
 % main equations for two variables:  L1, Fm
[Kf, alpha, beta, ks, tau, Vp0, Vd, F0, Fp0, phi, gamma, rcs, rfs, d0, Npod]...
    = deal(spara(1),spara(2), spara(3),spara(4),spara(5),spara(6),spara(7)...
    ,spara(8),spara(9), spara(10),spara(11),spara(12),spara(13), spara(14), spara(15));

dydt = zeros(3*Npod,1);
randnoiseM = randn(Npod,1);
randnoiseP = randn(Npod,1);
NoiseVarM = 0*100;
NoiseVarP = 0*Vp0*0.01;
mu = 0.06;
Dc = d0^2/80;
%need rearrange the equations below, as dydt is also involved in dFmdt
for i=1:Npod
    dFmdt = (F0+(alpha-gamma/Kf/rfs)*ks*y(3*i-2)/phi-(1-gamma/Kf)*y(3*i-1))/tau...
        +NoiseVarM*randnoiseM(i)/tau;
    % recalculate the protrusive force based on the podosome connectivity
    neighbors = podoconnect{i};
    Fp0n = Fp0;
    Vp0n = Vp0+y(3*i)*beta;
    dcdx = 0;

    if length(neighbors)==6
       for j=1:3
         jn1 = neighbors(j);
         jn2 = neighbors(7-j);
         dcdx = dcdx+(y(3*jn1)+y(3*jn2)-2*y(3*i))/d0^2;
       
       end
    end
    % the three main governing equation for each podosome i
    dydt(3*i-2) = (Vp0n*(1-ks*y(3*i-2)/Fp0n/rcs)-Vd+1/Kf*phi*dFmdt...
        +NoiseVarP*randnoiseP(i))/(ks/Kf/rfs+1);
    dydt(3*i-1) = (F0+(alpha-gamma/Kf/rfs)*ks*y(3*i-2)/phi-(1-gamma/Kf)*y(3*i-1))/tau...
        +NoiseVarM*randnoiseM(i)/tau;
    dydt(3*i) = Dc*dcdx-mu*(Vp0n*(1-ks*y(3*i-2)/Fp0n/rcs)-Vd);
end

end

