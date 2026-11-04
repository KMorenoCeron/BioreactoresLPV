%  BLOQUE 2 - OPTIMIZADOR POR COMBINACION OBSERVABLE w  (validacion)

%  Idea clave: el optimo NO depende de Sin y Ain por separado, solo de
%  su combinacion  w = Ain + k*Sin, que SI es observable desde A:
%
%        A = Ain + k*(Sin - S*)   =>   w = Ain + k*Sin = A + k*S*
%
%  Con w se calcula todo el punto de operacion optimo, sin necesitar
%  conocer Ain (que es una perturbacion desconocida):
%
%        D_op = (mu_max/alpha)*(1 - sqrt(k*kS/(w + k*kS)))
%        A*   = w - k*S*(D_op)
%        Z*   = (Sin_hat - S*)/alpha   (unica que usa un reparto; solo
%                                       afecta el transitorio del control)
%
%  LIMITE DE WASHOUT: la biomasa se lava si D supera mu_max/alpha=1.9.
%  Se clampa w a su envelope fisico (mata picos transitorios de w_hat)
%  y D_op a [D_min, D_max], con D_max < 1.9 por margen, y D_min como
%  piso de excitacion del observador.
%
%  Montaje de la prueba: Ain real FIJO fuera de nominal (=35, desconocido)
%  y Sin con escalones. Se valida que w_hat sigue a w real y que la
%  productividad P_A alcanza P_A,max, con D_op siempre lejos del washout.
%  Sin funciones. Todo inline. Incluye el observador del Bloque 1.



clear; clc;

%%  Parametros 
alpha  = 0.16;  mu_max = 0.38;  kS = 24.0;  k = 3.80;
Ain_nom= 50.0;                 % Ain NOMINAL que asume el observador
tau    = 5.76/(60*24);
p1 = -1.5;  p2 = -2.5;  l2 = p1*p2/k;      % l1 depende de D (se calcula en el lazo)
Dwo    = mu_max/alpha;          % dilucion de washout (1.9)

%%  Clamps (limite de washout y excitacion) 
w_min = 90.0;   w_max = 460.0;  % envelope fisico de w = Ain + k*Sin
D_min = 0.10;   D_max = 1.20;   % piso (excitacion)  y  tope (< washout 1.9)

%%  Perturbaciones de la prueba 
Ain_real = 35.0;                % Ain real, FIJO y DESCONOCIDO (!= nominal)
% perfil de Sin: escalones  45 -> 75 -> 55

%%  Condiciones iniciales (equilibrio del primer tramo) 
Sin0 = 60;
Dt0  = (mu_max/alpha)*(1 - sqrt(k*kS/(Ain_real + k*(kS+Sin0))));
S0   = alpha*kS*Dt0/(mu_max - alpha*Dt0);
Z0   = (Sin0 - S0)/alpha;
A0   = Ain_real + k*(Sin0 - S0);
mu0  = mu_max*S0/(kS+S0);

%%  Simulacion 
Tfin = 180;  N = round(Tfin/tau);  t = (0:N-1)*tau;
Z=Z0; S=S0; A=A0;                         % planta real
A_hat=0.5*A0; phi_hat=0.5*mu0*Z0; Z_hat=0.5*Z0;   % observador (50% error)
% D_op inicial con el clamp aplicado
w_ini = A0 + k*S0;  w_ini = min(max(w_ini,w_min),w_max);
D = (mu_max/alpha)*(1 - sqrt(k*kS/(w_ini + k*kS)));  D = min(max(D,D_min),D_max);

LOG=struct('Sin',zeros(1,N),'Sinh',zeros(1,N),'wr',zeros(1,N),'wh',zeros(1,N), ...
           'D',zeros(1,N),'PA',zeros(1,N),'PAmax',zeros(1,N));

for i = 1:N
    tk = t(i);
    if     tk<60,  Sin=45;
    elseif tk<120, Sin=75;
    else           Sin=55;  end

    %  PLANTA real 
    mu = mu_max*S/(kS+S);
    Z = Z + tau*((mu - alpha*D)*Z);
    S = S + tau*((Sin - S)*D - mu*Z);
    A = A + tau*((Ain_real - A)*D + k*mu*Z);

    %  BLOQUE 1: observador PI 
    l1 = -(p1+p2) - D;                        % gana con el D actual
    e  = A - A_hat;
    A_hat   = A_hat   + tau*((Ain_nom - A_hat)*D + k*phi_hat + l1*e);
    phi_hat = phi_hat + tau*(l2*e);
    Z_hat   = Z_hat   + tau*(phi_hat - alpha*D*Z_hat);
    mu_hat  = min(phi_hat/max(Z_hat,1e-6), 0.999*mu_max);
    S_hat   = min(max(kS*mu_hat/(mu_max - mu_hat),0),100);

    %  BLOQUE 2: combinacion observable y optimo, utiliza S_hat
    w_hat = A + k*S_hat;                      % w = Ain + k*Sin, observable
    w_hat = min(max(w_hat, w_min), w_max);    % clamp fisico (mata picos)
    D_op  = (mu_max/alpha)*(1 - sqrt(k*kS/(w_hat + k*kS)));
    D_op  = min(max(D_op, D_min), D_max);     % clamp washout + excitacion
    S_op  = alpha*kS*D_op/(mu_max - alpha*D_op);
    A_op  = w_hat - k*S_op;                   % referencia de acetato (solo w)
    % Z_op usaria Sin_hat (reparto), aqui no hace falta para validar el optimo
    D = D_op;                                 % feedforward (sin controlador aun)

    %  Referencias reales para comparar 
    w_real = Ain_real + k*Sin;
    Dtrue  = (mu_max/alpha)*(1 - sqrt(k*kS/(w_real + k*kS)));
    Strue  = alpha*kS*Dtrue/(mu_max - alpha*Dtrue);
    Atrue  = Ain_real + k*(Sin - Strue);

    LOG.Sin(i)=Sin;  LOG.Sinh(i)=S_hat + phi_hat/D;  % Sin_hat (sesgado, informativo)
    LOG.wr(i)=w_real; LOG.wh(i)=w_hat;
    LOG.D(i)=D_op;   LOG.PA(i)=A*D_op;  LOG.PAmax(i)=Atrue*Dtrue;
end

%%  Validacion 
m = t>30;   % tras convergencia inicial del observador
fprintf('Ain real=%.0f (desconocido).  Tras t>30 d:\n', Ain_real);
fprintf('  error medio de w  = %.3f %%\n', 100*mean(abs(LOG.wh(m)-LOG.wr(m))./LOG.wr(m)));
fprintf('  eficiencia media  = %.2f %%  (P_A / P_A,max)\n', 100*mean(LOG.PA(m)./LOG.PAmax(m)));
fprintf('  D_op maximo visto = %.4f   (washout = %.2f)\n', max(LOG.D), Dwo);

%% ---------- Graficas ----------
figure('Color','w','Name','Bloque 2 - optimizador por w');
subplot(4,1,1); plot(t,LOG.Sin,'k-',t,LOG.Sinh,'b--','LineWidth',1.2); grid on;
ylabel('S_{in}'); legend('S_{in} real','S_{in} est (sesgado)','Location','best');
title('Bloque 2: A_{in} real=35 desconocido, S_{in} con escalones');
subplot(4,1,2); plot(t,LOG.wr,'k-',t,LOG.wh,'r--','LineWidth',1.2); grid on;
ylabel('w'); legend('w real','w_{hat}=A+k S_{hat}','Location','best');
subplot(4,1,3); plot(t,LOG.D,'g-','LineWidth',1.3); hold on;
yline(Dwo,'r:','washout 1.9'); yline(D_max,'--','tope 1.2'); grid on;
ylabel('D_{op}'); ylim([0 2]); legend('D_{op} (con clamp)','Location','best');
subplot(4,1,4); plot(t,LOG.PAmax,'k-',t,LOG.PA,'g--','LineWidth',1.2); grid on;
ylabel('P_A'); xlabel('t [d]'); legend('P_{A,max}','P_A logrado','Location','best');
