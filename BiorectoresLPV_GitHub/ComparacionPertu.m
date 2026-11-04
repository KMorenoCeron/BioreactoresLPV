
%  COMPARACION DE PERTURBACIONES:  varia Sin   vs   varia Ain
%  Base: lazo LPV enfoque w (alpha=0.16, ganancias afines).
%  El controlador SIEMPRE usa Ain_nom=50 (observador y Zop).
%   - Caso SIN: Ain real = 50 fijo, Sin sigue un perfil.
%   - Caso AIN: Sin = 60 fijo, Ain real sigue un perfil (mismos % de salto).

clear; clc; close all;

%%  Parametros planta
P.alpha=0.16; P.mu_max=0.38; P.kS=24.0; P.k=3.80;
P.Ain_nom=50.0;          % lo que CREE el controlador (fijo siempre)
Sin_nom =60.0;           % base del sustrato
tau=5.76/(60*24);
P.p1=-1.5; P.p2=-2.5;

%%  Clamps / estabilizadores
P.w_min=90; P.w_max=460; P.D_min=0.10; P.D_max=1.20; P.Dsat=[0.01 1.2];
P.t_ctrl=25; P.tau_ref=3.0; P.dD_up=0.8; P.dD_dn=0.05;

%%  Ganancias LPV afines
P.K0=[-0.000891035  3.72682e-09  0.0160904 -0.00703808];
P.K1=[ 0.0280023  8.03147e-13  5.98639e-10  2.32771e-11];
P.K2=[ 9.73774e-13  0.0280023  2.17105e-11  2.63466e-12];
P.R1=[0.0020 0.2400]; P.R2=[0.004 7.8236];

%%  Perfil de saltos RELATIVOS (mismos para Sin y Ain, mismos tiempos)
% rel = multiplicador sobre el valor base en cada tramo
tsw = [50 95 180 270];                 % tiempos de cambio
rel = [0.9 1.2 1.00 1.3 0.8];  % = 65/60,75/60,60/60,90/60,40/60
% -> Sin: 65 75 60 90 40   |   Ain: 54.2 62.5 50 75 33.3 (mismos %)

Tfin=365; N=round(Tfin/tau); t=(0:N-1)*tau;

% perfiles como funciones anonimas
tramo = @(tk) 1 + sum(tk>=tsw);        % indice de tramo 1..5
Sin_var = @(tk) Sin_nom   * rel(tramo(tk));   % varia Sin
Ain_var = @(tk) P.Ain_nom * rel(tramo(tk));   % varia Ain (mismos %)
cte     = @(v) @(tk) v;                        % perfil constante

%%  CASO A: varia Sin  (Ain real = 50)
RA = simular(t,tau,P, Sin_var, cte(P.Ain_nom), Sin_nom);
%%  CASO B: varia Ain  (Sin = 60)
RB = simular(t,tau,P, cte(Sin_nom), Ain_var, Sin_nom);

%%  TABLA
fprintf('\n============ COMPARACION DE PERTURBACIONES ============\n');
fprintf('%-28s | %10s | %10s\n','Metrica (t>=25 d)','VARIA Sin','VARIA Ain');
fprintf('%-28s | %10.2f | %10.2f\n','Eficiencia P_A/P_Amax [%]',RA.eff,RB.eff);
fprintf('%-28s | %10.3f | %10.3f\n','TV(D)',RA.tvd,RB.tvd);
fprintf('%-28s | %10.3f | %10.3f\n','RMS(A - A*)',RA.rmsA,RB.rmsA);
fprintf('%-28s | %10.3f | %10.3f\n','Max |A - A*|',RA.maxA,RB.maxA);
fprintf('=======================================================\n');
if RA.eff>=RB.eff
  fprintf('>> Ain afecta MAS (menor eficiencia): caida de %.2f pts.\n',RA.eff-RB.eff);
else
  fprintf('>> Sin afecta MAS (menor eficiencia): caida de %.2f pts.\n',RB.eff-RA.eff);
end

%%  FIGURAS
figure('Color','w','Name','Comparacion Sin vs Ain');
subplot(3,1,1);
plot(t,RA.D,'b-','LineWidth',1.2); hold on; plot(t,RB.D,'r-','LineWidth',1.2);
yline(1.2,':k'); grid on; ylabel('D (d^{-1})'); ylim([0.4 1.4]);
legend('varia Sin','varia Ain','Location','best'); title('Accion de control D');
subplot(3,1,2);
plot(t,RA.PA./RA.PAm,'b-','LineWidth',1.2); hold on;
plot(t,RB.PA./RB.PAm,'r-','LineWidth',1.2); grid on;
ylabel('P_A / P_{A,max}'); ylim([0 3]); title('Eficiencia instantanea');
subplot(3,1,3);
plot(t,RA.A-RA.Aop,'b-','LineWidth',1.0); hold on;
plot(t,RB.A-RB.Aop,'r-','LineWidth',1.0); yline(0,':k'); grid on;
ylabel('A - A^*'); xlabel('t [d]'); title('Error de acetato vs objetivo');

%% =====================================================================
function R = simular(t,tau,P, Sin_fun, Ain_fun, Sin0)
  a=P.alpha; mm=P.mu_max; kS=P.kS; k=P.k; Ain_nom=P.Ain_nom;
  N=numel(t); l2=P.p1*P.p2/k;

  % CI con el valor base y Ain nominal
  D0=(mm/a)*(1-sqrt(k*kS/(Ain_nom+k*(kS+Sin0))));
  S0=a*kS*D0/(mm-a*D0); Z0=(Sin0-S0)/a; A0=Ain_nom+k*(Sin0-S0);
  mu0=mm*S0/(kS+S0);

  Z=Z0; S=S0; A=A0;
  A_hat=A0; phi_hat=0.5*mu0*Z0; Z_hat=0.5*Z0;
  w_ini=min(max(A0+k*S0,P.w_min),P.w_max);
  D=min(max((mm/a)*(1-sqrt(k*kS/(w_ini+k*kS))),P.D_min),P.D_max);
  Dopf=D0; Sopf=S0; Aopf=A0; Zopf=Z0; zI=0;

  R.D=zeros(1,N); R.A=zeros(1,N); R.Aop=zeros(1,N);
  R.PA=zeros(1,N); R.PAm=zeros(1,N);

  for i=1:N
    tk=t(i);
    Sin     = Sin_fun(tk);     % entra SOLO a la planta
    Ain_tru = Ain_fun(tk);     % Ain REAL (el controlador no lo sabe)

    % PLANTA (Euler) con Ain real
    mu=mm*S/(kS+S);
    Z=Z+tau*((mu-a*D)*Z);
    S=S+tau*((Sin-S)*D-mu*Z);
    A=A+tau*((Ain_tru-A)*D+k*mu*Z);

    % OBSERVADOR PI  (usa Ain_nom, no el real)
    l1=-(P.p1+P.p2)-D; e=A-A_hat;
    A_hat=A_hat+tau*((Ain_nom-A_hat)*D+k*phi_hat+l1*e);
    phi_hat=phi_hat+tau*(l2*e);
    Z_hat=Z_hat+tau*(phi_hat-a*D*Z_hat);
    mu_hat=min(phi_hat/max(Z_hat,1e-6),0.999*mm);
    S_hat=min(max(kS*mu_hat/(mm-mu_hat),0),100);

    % OPTIMO por w  (Zop usa Ain_nom)
    w=min(max(A+k*S_hat,P.w_min),P.w_max);
    Dop=min(max((mm/a)*(1-sqrt(k*kS/(w+k*kS))),P.D_min),P.D_max);
    Sop=a*kS*Dop/(mm-a*Dop);
    Aop=w-k*Sop;
    Zop=max(((w-Ain_nom)/k-Sop)/a,0);

    % scheduling
    rho1=min(max(mm*S_hat/(kS+S_hat),P.R1(1)),P.R1(2));
    rho2=min(max(mm*kS/(kS+S_hat)^2*Z_hat,P.R2(1)),P.R2(2));

    % gobernador
    c=tau/P.tau_ref;
    Dopf=Dopf+c*(Dop-Dopf); Sopf=Sopf+c*(Sop-Sopf);
    Aopf=Aopf+c*(Aop-Aopf); Zopf=Zopf+c*(Zop-Zopf);

    % controlador LPV + anti-windup
    D_prev=D;
    if tk>=P.t_ctrl
      Krho=P.K0+P.K1*rho1+P.K2*rho2;
      xaug=[Z_hat-Zopf; S_hat-Sopf; A-Aopf; zI];
      D=min(max(Dopf+Krho*xaug,P.Dsat(1)),P.Dsat(2));
      eA=A-Aopf;
      ns=(D>P.Dsat(1))&&(D<P.Dsat(2));
      sp=(D<=P.Dsat(1))&&(eA>0);
      st=(D>=P.Dsat(2))&&(eA<0);
      if ns||sp||st, zI=zI-eA*tau; end
    else
      D=Dopf;
    end
    D=min(max(D,D_prev-P.dD_dn*tau),D_prev+P.dD_up*tau);

    % REFERENCIA VERDADERA  (usa Ain real y Sin real)
    Dtr=(mm/a)*(1-sqrt(k*kS/(Ain_tru+k*(kS+Sin))));
    Str=a*kS*Dtr/(mm-a*Dtr);
    Atr=Ain_tru+k*(Sin-Str);

    R.D(i)=D; R.A(i)=A; R.Aop(i)=Aopf;
    R.PA(i)=A*D; R.PAm(i)=Atr*Dtr;
  end

  m=t>=P.t_ctrl;
  R.eff =100*mean(R.PA(m)./R.PAm(m));
  R.tvd =sum(abs(diff(R.D(m))));
  R.rmsA=sqrt(mean((R.A(m)-R.Aop(m)).^2));
  R.maxA=max(abs(R.A(m)-R.Aop(m)));
end