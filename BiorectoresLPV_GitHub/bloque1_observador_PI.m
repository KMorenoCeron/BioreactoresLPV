
%  BLOQUE 1 - OBSERVADOR PI DE TASA  (validacion aislada)
clear; clc;

%% Parametros del modelo 
alpha  = 0.16;      % fraccion de biomasa en fase liquida
mu_max = 0.38;      % tasa maxima de crecimiento [1/d]
kS     = 24.0;      % constante de semisaturacion [g/L]
k      = 3.80;      % coef. de rendimiento de acetato
Ain    = 50.0;      % acetato de entrada (nominal) [mmol/L]
tau    = 5.76/(60*24);   % paso de muestreo [d]  (5.76 min)

%% Punto de operacion OPTIMO del paper 
Sin  = 60.0;        % sustrato de entrada nominal [g/L]
Dop  = (mu_max/alpha)*(1 - sqrt(k*kS/(Ain + k*(kS+Sin))));  % dilucion optima
S_eq = alpha*kS*Dop/(mu_max - alpha*Dop);     % sustrato en el optimo
mu_eq= mu_max*S_eq/(kS+S_eq);                 % tasa de crecimiento en el optimo
Z_eq = (Sin - S_eq)/alpha;                    % biomasa en el optimo
A_eq = Ain + k*(Sin - S_eq);                  % acetato en el optimo
fprintf('Punto optimo:  Dop=%.4f  S*=%.3f  Z*=%.3f  A*=%.3f\n', Dop,S_eq,Z_eq,A_eq);

%% Ganancias del observador PI (ubicacion de polos) 
%  Se ubican los polos del error en p1, p2:
p1 = -1.5;  p2 = -2.5;
l1 = -(p1+p2) - Dop;      % gana la correccion proporcional sobre A_hat
l2 =  p1*p2 / k;          % gana la parte integral que estima phi
fprintf('Ganancias PI:  l1=%.4f  l2=%.4f\n\n', l1, l2);

%% Configuracion de la simulacion
Tfin = 60;  N = round(Tfin/tau);  t = (0:N-1)*tau;
D    = Dop;                 % dilucion fija en el optimo (lazo abierto)
Smax = 100.0;               % S fisico no supera el Sin maximo de operacion

% Estado REAL de la planta: arranca en el equilibrio
Z = Z_eq;  S = S_eq;  A = A_eq;
% Estado del OBSERVADOR: arranca con 50% de error
A_hat   = 0.75*A_eq;
phi_hat = 0.75*mu_eq*Z_eq;
Z_hat   = 0.75*Z_eq;

% Buffers para graficar
LOG_Z=zeros(1,N); LOG_Zh=zeros(1,N);
LOG_S=zeros(1,N); LOG_Sh=zeros(1,N);
LOG_A=zeros(1,N); LOG_Ah=zeros(1,N);

%%  Lazo 
for i = 1:N
    %  PLANTA real (un paso de Euler) 
    mu = mu_max*S/(kS+S);
    Z = Z + tau*((mu - alpha*D)*Z);
    S = S + tau*((Sin - S)*D - mu*Z);
    A = A + tau*((Ain - A)*D + k*mu*Z);

    %  OBSERVADOR PI (un paso) 
    e       = A - A_hat;                                   % error de acetato
    A_hat   = A_hat   + tau*((Ain - A_hat)*D + k*phi_hat + l1*e);  % P: corrige A_hat
    phi_hat = phi_hat + tau*(l2*e);                        % I: estima la tasa phi
    Z_hat   = Z_hat   + tau*(phi_hat - alpha*D*Z_hat);     % biomasa desde la tasa
    mu_hat  = min(phi_hat/max(Z_hat,1e-6), 0.999*mu_max);  % tasa de crecimiento
    S_hat   = kS*mu_hat/(mu_max - mu_hat);                 % sustrato (inversion Monod)
    S_hat   = min(max(S_hat,0.0), Smax);                   % clamp fisico (evita el pico)

    %  Registro 
    LOG_Z(i)=Z;  LOG_Zh(i)=Z_hat;
    LOG_S(i)=S;  LOG_Sh(i)=S_hat;
    LOG_A(i)=A;  LOG_Ah(i)=A_hat;
end

%%  Validacion: tiempo de convergencia (banda 2%) 
idxZ = find(abs(LOG_Zh-LOG_Z)./abs(LOG_Z) > 0.02, 1, 'last');
idxS = find(abs(LOG_Sh-LOG_S)./abs(LOG_S) > 0.02, 1, 'last');
tZ = (isempty(idxZ))*0 + (~isempty(idxZ))*t(max(idxZ,1));
tS = (isempty(idxS))*0 + (~isempty(idxS))*t(max(idxS,1));
fprintf('Convergencia (banda 2%%):  Z en %.1f d   S en %.1f d\n', tZ, tS);
fprintf('Error final:  Z = %.3f%%   S = %.3f%%\n', ...
        100*abs(LOG_Zh(end)-LOG_Z(end))/LOG_Z(end), ...
        100*abs(LOG_Sh(end)-LOG_S(end))/LOG_S(end));
%%  Graficas 
figure('Color','w','Name','Bloque 1 - Observador PI');
subplot(3,1,1);
plot(t,LOG_Z,'r--','LineWidth',1.3); hold on; plot(t,LOG_Zh,'b-','LineWidth',1.3);
grid on; ylabel('Z'); legend('Z real','Z est','Location','best');
title('Bloque 1 - Observador PI: convergencia desde 50% de error');
subplot(3,1,2);
plot(t,LOG_S,'r--','LineWidth',1.3); hold on; plot(t,LOG_Sh,'b-','LineWidth',1.3);
grid on; ylabel('S'); legend('S real','S est','Location','best');
subplot(3,1,3);
plot(t,LOG_A,'r--','LineWidth',1.3); hold on; plot(t,LOG_Ah,'b-','LineWidth',1.3);
grid on; ylabel('A'); xlabel('t [d]'); legend('A real','A est','Location','best');
