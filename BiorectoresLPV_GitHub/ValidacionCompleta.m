
%  LAZO CERRADO LPV - ENFOQUE w  (modular, sin funciones)
%  Observador PI -> optimizador por w -> scheduling -> controlador LPV
%   El punto de operacion sale de  w = A + k*S_hat  (NO divide por D),
%   Perfil del paper (65/75/60/90/40), 365 d, Ain=50.

clear; clc;

%%  Parametros 
alpha=0.16; mu_max=0.38; kS=24.0; k=3.80; Ain=50.0;
tau=5.76/(60*24);
p1=-1.5; p2=-2.5; l2=p1*p2/k;          % l1 depende de D (se calcula en el lazo)

%%  Clamps del optimizador 
w_min=90; w_max=460;  D_min=0.10; D_max=1.20;
Dsat=[0.01 1.2];

%%  Estabilizadores 
t_ctrl=25;        % el LPV entra tras converger el observador
tau_ref=3.0;      % gobernador de referencia
dD_up=0.8;        % slew: D SUBE rapido
dD_dn=0.05;       % slew: D BAJA lento  (clave para no desplomarse)

%%  Ganancias LPV afines: K(rho)=K0+K1*rho1+K2*rho2 
if exist('lpv_gains.mat','file')==2
    G=load('lpv_gains.mat');
    K0=G.K0; K1=G.K1; K2=G.K2; R1=G.R1; R2=G.R2;
    fprintf('Ganancias cargadas de lpv_gains.mat\n');
else
    K0=[-0.000891035  3.72682e-09  0.0160904 -0.00703808];
    K1=[ 0.0280023  8.03147e-13  5.98639e-10  2.32771e-11];
    K2=[ 9.73774e-13  0.0280023  2.17105e-11  2.63466e-12];
    R1=[0.0020 0.2400]; R2=[0.004 7.8236];
    warning('lpv_gains.mat no encontrado: usando ganancias del ultimo diseno.');
end

%%  Condiciones iniciales 
Sin0=60;
D0=(mu_max/alpha)*(1-sqrt(k*kS/(Ain+k*(kS+Sin0))));
S0=alpha*kS*D0/(mu_max-alpha*D0); Z0=(Sin0-S0)/alpha; A0=Ain+k*(Sin0-S0);
mu0=mu_max*S0/(kS+S0);

Tfin=365; N=round(Tfin/tau); t=(0:N-1)*tau;
Z=Z0; S=S0; A=A0;                               % planta real
A_hat=A0; phi_hat=0.5*mu0*Z0; Z_hat=0.5*Z0;     % observador
w_ini=min(max(A0+k*S0,w_min),w_max);
D=min(max((mu_max/alpha)*(1-sqrt(k*kS/(w_ini+k*kS))),D_min),D_max);
Dopf=D0; Sopf=S0; Aopf=A0; Zopf=Z0; zI=0.0;

LA=zeros(1,N);LAop=zeros(1,N);LD=zeros(1,N);LPA=zeros(1,N);LPAm=zeros(1,N);
LZ=zeros(1,N);LZh=zeros(1,N);LS=zeros(1,N);LSh=zeros(1,N);Lr1=zeros(1,N);Lr2=zeros(1,N);
LDop=zeros(1,N);
%%  Lazo 
for i=1:N
    tk=t(i);
    if     tk<50,  Sin=60;
    elseif tk<95,  Sin=55;
    elseif tk<180, Sin=60;
    elseif tk<270, Sin=80;
    else           Sin=90;  end

    %  PLANTA (Euler) 
    mu=mu_max*S/(kS+S);
    Z=Z+tau*((mu-alpha*D)*Z);
    S=S+tau*((Sin-S)*D-mu*Z);
    A=A+tau*((Ain-A)*D+k*mu*Z);

    %  BLOQUE 1: observador PI 
    l1=-(p1+p2)-D;  e=A-A_hat;
    A_hat=A_hat+tau*((Ain-A_hat)*D+k*phi_hat+l1*e);
    phi_hat=phi_hat+tau*(l2*e);
    Z_hat=Z_hat+tau*(phi_hat-alpha*D*Z_hat);
    mu_hat=min(phi_hat/max(Z_hat,1e-6),0.999*mu_max);
    S_hat=min(max(kS*mu_hat/(mu_max-mu_hat),0),100);

    %  BLOQUE 2: optimo por w (sin phi/D) 
    w=min(max(A+k*S_hat,w_min),w_max);
    Dop=min(max((mu_max/alpha)*(1-sqrt(k*kS/(w+k*kS))),D_min),D_max);
    Sop=alpha*kS*Dop/(mu_max-alpha*Dop);
    Aop=w-k*Sop;
    Zop=max(((w-Ain)/k - Sop)/alpha, 0);        % Z* robusto (reparto con Ain nominal)

    %  BLOQUE 3b: scheduling 
    rho1=min(max(mu_max*S_hat/(kS+S_hat),R1(1)),R1(2));
    rho2=min(max(mu_max*kS/(kS+S_hat)^2*Z_hat,R2(1)),R2(2));

    %  Gobernador de referencia 
    c=tau/tau_ref;
    Dopf=Dopf+c*(Dop-Dopf); Sopf=Sopf+c*(Sop-Sopf);
    Aopf=Aopf+c*(Aop-Aopf); Zopf=Zopf+c*(Zop-Zopf);

    %  BLOQUE 4: controlador LPV con integrador + anti-windup 
    D_prev=D;
    if tk>=t_ctrl
        Krho=K0 + K1*rho1 + K2*rho2;            % ganancia afin (forma PLC)
        xaug=[Z_hat-Zopf; S_hat-Sopf; A-Aopf; zI];
        D=min(max(Dopf+Krho*xaug,Dsat(1)),Dsat(2));
        eA=A-Aopf;
        no_sat   =(D>Dsat(1))&&(D<Dsat(2));
        sale_piso=(D<=Dsat(1))&&(eA>0);
        sale_tope=(D>=Dsat(2))&&(eA<0);
        if no_sat||sale_piso||sale_tope, zI=zI-eA*tau; end
    else
        D=Dopf;
    end
    % slew ASIMETRICO: sube rapido, baja lento
    D=min(max(D, D_prev-dD_dn*tau), D_prev+dD_up*tau);

    %  referencia real para comparar 
    Dtrue=(mu_max/alpha)*(1-sqrt(k*kS/(Ain+k*(kS+Sin))));
    Strue=alpha*kS*Dtrue/(mu_max-alpha*Dtrue);
    Atrue=Ain+k*(Sin-Strue);

    LA(i)=A; LAop(i)=Aopf; LD(i)=D; LPA(i)=A*D; LPAm(i)=Atrue*Dtrue;
    LZ(i)=Z; LZh(i)=Z_hat; LS(i)=S; LSh(i)=S_hat; Lr1(i)=rho1; Lr2(i)=rho2;
    LDop(i)=Dopf;
end

%%  Metricas 
m=t>=t_ctrl;
fprintf('Eficiencia media P_A/P_A,max (t>=%d) = %.2f %%\n', t_ctrl, 100*mean(LPA(m)./LPAm(m)));
fprintf('D en [%.3f, %.3f]   (washout = %.2f)\n', min(LD), max(LD), mu_max/alpha);
fprintf('rho1 in [%.4f, %.4f]   rho2 in [%.4f, %.4f]\n', min(Lr1),max(Lr1),min(Lr2),max(Lr2));
fprintf('Variacion total TV(D) = %.2f\n', sum(abs(diff(LD(m)))));

%%  Graficas 
figure('Color','w','Name','Lazo cerrado LPV - enfoque w');
subplot(4,1,1); plot(t,LD,'b-','LineWidth',1.2); hold on; yline(1.2,':k','u_{max}'); grid on;
ylabel('D (d^{-1})'); ylim([0.4 1.4]); title('Tasa de dilucion');
subplot(4,1,2); plot(t,LPAm,'k-','LineWidth',1.0); hold on; plot(t,LPA,'b-','LineWidth',1.2); grid on;
ylabel('P_A'); ylim([0 350]); legend('P_{A,max}','P_A','Location','best'); title('Productividad');
subplot(4,1,3); plot(t,LA,'b-','LineWidth',1.2); hold on; plot(t,LAop,'r--','LineWidth',1.0); grid on;
ylabel('A'); legend('A','A^*','Location','best'); title('Acetato vs objetivo');
subplot(4,1,4); plot(t,LZ,'r--',t,LZh,'b-',t,LS,'m--',t,LSh,'g-','LineWidth',1.0); grid on;
ylabel('Z, S'); xlabel('t [d]'); legend('Z','Z_{hat}','S','S_{hat}','Location','best'); title('Estados real vs estimado');

%%
figure('Color','w','Name','Accion de control D');
subplot(2,1,1);
plot(t,LD,'b-','LineWidth',1.3); hold on;
plot(t,LDop,'r--','LineWidth',1.0);
yline(1.2,':k','u_{max}'); grid on;
ylabel('D (d^{-1})'); ylim([0.4 1.4]);
legend('D aplicada','D_{op} (feedforward)','Location','best');
title('Accion de control D');
subplot(2,1,2);
plot(t,LD-LDop,'m-','LineWidth',1.0); yline(0,':k'); grid on;
ylabel('K(\rho)\cdot x_{aug}'); xlabel('t [d]');
title('Correccion del controlador (feedback) = D - D_{op}');