
%  GRAFICAS DE COMPARACION (LPV, enfoque w) 
%  Lazo cerrado con el optimizador por w (A + k*S_hat, sin phi/D),
%     FIG 1 : ESTIMACION   -> Z, S, A  real vs estimado/objetivo
%     FIG 2 : PRODUCCION   -> P_A vs P_A,max
%     FIG 3 : ACCION DE CONTROL -> D aplicada, D_op y correccion
%  Perfil de Sin (~365 d): 65/75/60/90/50.

clear; clc; close all;

%%  Parametros 
par.alpha=0.16; par.mu_max=0.38; par.kS=24.0; par.k=3.80; par.Ain=50.0;
tau = 5.76/(60*24);

%%  Ganancias PI de tasa (Bloque 1) 
Sin_nom=60.0;
[Dop_nom,~,~,~] = punto_operacion_true(Sin_nom,par);
p_est=[-1.5 -2.5];
obs.l1=-(p_est(1)+p_est(2))-Dop_nom;  obs.l2=(p_est(1)*p_est(2))/par.k;

%%  Ganancias LPV (Bloque 3) 
if exist('lpv_gains.mat','file')==2
    S = load('lpv_gains.mat');
    K.K11=S.K11; K.K21=S.K21; K.K12=S.K12; K.K22=S.K22;
    R1 = S.R1;   R2 = S.R2;
    fprintf('Ganancias LPV cargadas de lpv_gains.mat.\n');
else
    K.K11 = [-0.0028591  0.00098845  0.012187  -0.009504];
    K.K21 = [ 0.0038034  0.00098539  0.012196  -0.0095106];
    K.K12 = [-0.0028598  0.2187      0.012181  -0.0095049];
    K.K22 = [ 0.0038028  0.21871     0.012199  -0.0095123];
    R1=[0.0020 0.2400];   R2=[0.0485 7.8236];
    warning('lpv_gains.mat no encontrado: usando ganancias historicas.');
end
Dsat=[0.01 1.2];

%%  Clamps del optimizador y slew, permite que no bajé de golpe, representativo del sistema
w_min=90; w_max=460;  D_min=0.10; D_max=1.20;
t_ctrl=25.0; tau_ref=3.0;
dD_up=0.8; dD_dn=0.05;        % D sube rapido, baja lento (clave)

%%  Config simulacion 
Tfin=365; N=round(Tfin/tau); t=(0:N-1)*tau;

% Condiciones iniciales
Sin0=perfil_Sin(0); [D0,S0,A0,Z0]=punto_operacion_true(Sin0,par);
x=[Z0;S0;A0];
piA=A0; piPhi=0.5*mu_op_of(S0,par)*Z0; piZ=0.5*Z0;
Dopf=D0; Sopf=S0; Aopf=A0; Zopf=Z0;
D=D0; zI=0.0;

% Buffers
Sin_=zeros(1,N); A_=zeros(1,N); D_=zeros(1,N); Dop_=zeros(1,N);
S_=zeros(1,N); Z_=zeros(1,N); Zh_=zeros(1,N); Sh_=zeros(1,N);
Aop_=zeros(1,N); PA_=zeros(1,N); PAmax_=zeros(1,N);

%%  Lazo 
for kk=1:N
    tk=t(kk); Sin_true=perfil_Sin(tk);
    x = x + tau*modelo_df(x,D,Sin_true,par);  A=x(3);

    % Bloque 1: observador PI
    [piA,piPhi,piZ,S_hat,~]=observador_tasa_PI(piA,piPhi,piZ,D,A,par,obs,tau);
    Z_hat=piZ;

    % Bloque 3b: scheduling
    [rho1,rho2]=parametros_lpv(S_hat,Z_hat,par,R1,R2);

    % Bloque 2: optimo por w (sin phi/D)
    [Dop,Sop,Aop,Zop]=punto_operacion_w(A,S_hat,par,w_min,w_max,D_min,D_max);

    % Gobernador de referencia
    c=tau/tau_ref;
    Dopf=Dopf+c*(Dop-Dopf); Sopf=Sopf+c*(Sop-Sopf);
    Aopf=Aopf+c*(Aop-Aopf); Zopf=Zopf+c*(Zop-Zopf);

    % Bloque 4: controlador LPV
    D_prev=D;
    if tk>=t_ctrl
        [D,zI]=controlador_lpv(Dopf,Zopf,Sopf,Aopf,Z_hat,S_hat,A,rho1,rho2,zI,K,R1,R2,Dsat,tau);
    else
        D=Dopf;
    end
    % slew ASIMETRICO
    D=min(max(D,D_prev-dD_dn*tau),D_prev+dD_up*tau);

    % Log
    Sin_(kk)=Sin_true; A_(kk)=A; D_(kk)=D; Dop_(kk)=Dopf;
    S_(kk)=x(2); Z_(kk)=x(1); Zh_(kk)=Z_hat; Sh_(kk)=S_hat;
    Aop_(kk)=Aopf; PA_(kk)=A*D;
    [Dm,~,Am,~]=punto_operacion_true(Sin_true,par); PAmax_(kk)=Am*Dm;
end

%%  Metricas 
m=t>=t_ctrl;
fprintf('Eficiencia media P_A/P_A,max = %.2f %%\n', 100*mean(PA_(m)./PAmax_(m)));
fprintf('TV(D) = %.2f   D en [%.3f, %.3f]\n', sum(abs(diff(D_(m)))), min(D_), max(D_));

%%  Curvas del paper (DIGITALIZADAS a ojo, APROXIMADAS) 
tgrid=0:0.25:365;
tD_ST=[0 50 51 60 90 96 100 108 180 181 190 270 271 277 285 300 365];
vD_ST=[0.87 0.87 1.20 0.97 0.97 0.90 0.50 0.90 0.96 1.30 1.02 1.00 1.45 0.90 0.74 0.78 0.80];
D_ST=interp1(tD_ST,vD_ST,tgrid,'linear');
tP_ST=[0 10 40 50 56 70 92 96 105 175 182 200 265 272 285 300 365];
vP_ST=[0 55 160 190 250 235 228 150 172 173 260 300 305 150 92 110 110];
PA_ST=interp1(tP_ST,vP_ST,tgrid,'linear');


%% FIGURAS DE COMPARACION
azul=[0.12 0.44 0.70]; rojo=[0.85 0.33 0.31];

%  FIG 1: ESTIMACION 
figure('Color','w','Name','FIG1 Estimacion');
subplot(3,1,1);
plot(t,Z_,'--','Color',rojo,'LineWidth',1.2); hold on;
plot(t,Zh_,'-','Color',azul,'LineWidth',1.2); grid on; xlim([0 365]);
ylabel('Z'); legend('Z real','Z estimado','Location','best'); title('Estimacion de biomasa');
subplot(3,1,2);
plot(t,S_,'--','Color',rojo,'LineWidth',1.2); hold on;
plot(t,Sh_,'-','Color',azul,'LineWidth',1.2); grid on; xlim([0 365]); ylim([0 100]);
ylabel('S'); legend('S real','S estimado','Location','best'); title('Estimacion de sustrato (DQO)');
subplot(3,1,3);
plot(t,A_,'-','Color',azul,'LineWidth',1.2); hold on;
plot(t,Aop_,'--','Color',rojo,'LineWidth',1.0); grid on; xlim([0 365]);
ylabel('A'); xlabel('t [d]'); legend('A medido','A^* objetivo','Location','best'); title('Acetato vs objetivo');

%  FIG 2: PRODUCCION  (LPV vs optimo vs paper) 
figure('Color','w','Name','FIG2 Produccion');
plot(t,PAmax_,'k-','LineWidth',1.0); hold on;
plot(tgrid,PA_ST,'--','Color',rojo,'LineWidth',1.2);        % paper ST (digit.)
plot(t,PA_,'-','Color',azul,'LineWidth',1.3);               
grid on; xlim([0 365]); ylim([0 350]);
xlabel('t [d]'); ylabel('P_A (mmol L^{-1} d^{-1})');
legend('P_{A,max} (optimo)','P_{A,ST} (paper, digit.)','P_{A,LPV}','Location','best');
title('Productividad: LPV vs ST (paper) vs optimo');

% --- FIG 3: ACCION DE CONTROL ---
figure('Color','w','Name','FIG3 Accion de control');
subplot(2,1,1);
plot(t,D_,'-','Color',azul,'LineWidth',1.3); hold on;
plot(tgrid,D_ST,':','Color',rojo,'LineWidth',1.2);          % paper ST (digit.)
yline(1.2,':k','u_{max}'); grid on; xlim([0 365]); ylim([0.4 1.5]);
ylabel('D (d^{-1})'); legend('D_{LPV}','D_{ST} (paper, digit.)','Location','best');
title('Accion de control: LPV vs ST (paper)');
% ============================================================
%  FUNCIONES
% ============================================================
function dx=modelo_df(x,D,Sin,par)
    mu=par.mu_max*x(2)/(par.kS+x(2));
    dx=[-par.alpha*x(1)*D+mu*x(1); (Sin-x(2))*D-mu*x(1); (par.Ain-x(3))*D+par.k*mu*x(1)];
end

function [A_hat,phi_hat,Z_hat,S_hat,mu_hat]=observador_tasa_PI(A_hat,phi_hat,Z_hat,D,A_med,par,obs,tau)
    e=A_med-A_hat;
    A_hat=A_hat+tau*((par.Ain-A_hat)*D+par.k*phi_hat+obs.l1*e);
    phi_hat=phi_hat+tau*(obs.l2*e);
    Z_hat=Z_hat+tau*(phi_hat-par.alpha*D*Z_hat);
    mu_hat=min(phi_hat/max(Z_hat,1e-6),0.999*par.mu_max);
    S_hat=min(max(par.kS*mu_hat/(par.mu_max-mu_hat),0),100);
end

function [rho1,rho2]=parametros_lpv(S_hat,Z_hat,par,R1,R2)
    rho1=par.mu_max*S_hat/(par.kS+S_hat);
    rho2=par.mu_max*par.kS/(par.kS+S_hat)^2*Z_hat;
    rho1=min(max(rho1,R1(1)),R1(2)); rho2=min(max(rho2,R2(1)),R2(2));
end

function [D_op,S_op,A_op,Z_op]=punto_operacion_w(A,S_hat,par,w_min,w_max,D_min,D_max)
    % Optimo por la combinacion observable w = A + k*S_hat (NO divide por D)
    w = min(max(A + par.k*S_hat, w_min), w_max);
    D_op = (par.mu_max/par.alpha)*(1 - sqrt(par.k*par.kS/(w + par.k*par.kS)));
    D_op = min(max(D_op, D_min), D_max);
    S_op = par.alpha*par.kS*D_op/(par.mu_max - par.alpha*D_op);
    A_op = w - par.k*S_op;
    Z_op = max(((w - par.Ain)/par.k - S_op)/par.alpha, 0);   % reparto con Ain nominal
end

function [D,zI]=controlador_lpv(Dop,Zop,Sop,Aop,Zh,Sh,A,rho1,rho2,zI,K,R1,R2,Dsat,tau)
    a2=(rho1-R1(1))/(R1(2)-R1(1)); a1=1-a2;
    b2=(rho2-R2(1))/(R2(2)-R2(1)); b1=1-b2;
    Krho=a1*b1*K.K11 + a2*b1*K.K21 + a1*b2*K.K12 + a2*b2*K.K22;
    xaug=[Zh-Zop; Sh-Sop; A-Aop; zI];
    Dus=Dop + Krho*xaug;
    D=min(max(Dus,Dsat(1)),Dsat(2));
    eA=A-Aop;
    no_sat   = (D>Dsat(1)) && (D<Dsat(2));
    sale_piso= (D<=Dsat(1)) && (eA>0);
    sale_tope= (D>=Dsat(2)) && (eA<0);
    if no_sat || sale_piso || sale_tope
        zI=zI-eA*tau;
    end
end

function [D_op,S_op,A_op,Z_op]=punto_operacion_true(Sin,par)
    D_op=(par.mu_max/par.alpha)*(1-sqrt(par.k*par.kS/(par.Ain+par.k*(par.kS+Sin))));
    S_op=par.alpha*par.kS*D_op/(par.mu_max-par.alpha*D_op);
    A_op=par.Ain+par.k*Sin-par.k*S_op; Z_op=(Sin-S_op)/par.alpha;
end

function m=mu_op_of(S,par), m=par.mu_max*S/(par.kS+S); end

function Sin=perfil_Sin(t)
    if     t<50,   Sin=65;
    elseif t<95,   Sin=75;
    elseif t<180,  Sin=60;
    elseif t<270,  Sin=90;
    else,          Sin=50;
    end
end