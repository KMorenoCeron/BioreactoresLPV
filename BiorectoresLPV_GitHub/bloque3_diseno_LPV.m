clear; clc; yalmip('clear');
%% ====================================================================
%  DISENO CONTROLADOR LPV  -  H-infinito + Region D  (Bloques 3 y 4)
%  Biorreactor de fermentacion oscura.  Estructura afin / politopica.
%  --------------------------------------------------------------------
%  Mismo flujo que el codigo de referencia (sistema de estanques):
%   1. Jacobiano simbolico
%   2. Parametros y punto de operacion
%   3. Rangos de los parametros variantes (barrido) -> politopo
%   4. Modelo LPV afin: A(rho)=A0+rho1*A1+rho2*A2, B constante
%   5. Estructura ROLMIP (A, B, canal de perturbacion y de desempeno)
%   6. LMI: H-infinito (lema real acotado) + Region D (alpha, r, theta)
%   7. Ganancia K = Z*inv(P) y ganancias de los 4 vertices
%   8. Regresion afin: K(rho) = K0 + K1*rho1 + K2*rho2   (para el PLC)
%
%  OFFLINE: se corre una vez en PC. Entrega R1,R2 y K0,K1,K2 para
%  programar el gain scheduling en el PLC (ninguna LMI corre en el PLC).
%  Requiere YALMIP + ROLMIP + SeDuMi.
%% ====================================================================

%% 1) JACOBIANO SIMBOLICO (solo para mostrar)
syms sZ sS sA sD sSin sAin s_alpha s_mumax s_kS s_k real
sMu = s_mumax*sS/(s_kS+sS);
sf  = [(sMu - s_alpha*sD)*sZ;
       (sSin - sS)*sD - sMu*sZ;
       (sAin - sA)*sD + s_k*sMu*sZ];
Jx = jacobian(sf,[sZ;sS;sA]);
Jd = jacobian(sf,sD);
disp('Jacobiano df/dx:'); disp(Jx)
disp('Jacobiano df/dD:'); disp(Jd)
clear sZ sS sA sD sSin sAin s_alpha s_mumax s_kS s_k sMu sf Jx Jd

%% 2) PARAMETROS Y PUNTO DE OPERACION (optimo, Sin nominal)
alpha=0.16; mu_max=0.38; kS=24.0; k=3.80; Ain=50.0; Sin=60.0;
Dop = (mu_max/alpha)*(1 - sqrt(k*kS/(Ain + k*(kS+Sin))));
Sop = alpha*kS*Dop/(mu_max - alpha*Dop);
mu_op = mu_max*Sop/(kS+Sop);
Zop = (Sin - Sop)/alpha;
Aop = Ain + k*(Sin - Sop);
fprintf('Punto de operacion:  Dop=%.4f  S*=%.3f  Z*=%.3f  A*=%.3f\n', Dop,Sop,Zop,Aop);

%% 3) RANGOS DE rho1, rho2 (barrido de la region de operacion)
Dg   = linspace(0.01,1.2,60);
Sing = linspace(20,100,60);
r1=[]; r2=[];
for Dq = Dg
    Sx = alpha*kS*Dq/(mu_max-alpha*Dq);
    if Sx<=0, continue; end
    mux = mu_max*Sx/(kS+Sx);
    dmu = mu_max*kS/(kS+Sx)^2;
    for Sq = Sing
        if Sx>=Sq, continue; end
        Zx = (Sq - Sx)/alpha;
        r1(end+1) = mux;
        r2(end+1) = dmu*Zx;
    end
end
R1 = [min(r1) max(r1)];  R2 = [min(r2) max(r2)];
fprintf('rho1 in [%.4f, %.4f]\n', R1(1),R1(2));
fprintf('rho2 in [%.4f, %.4f]\n\n', R2(1),R2(2));
lim = [R1; R2];

%% 4) MODELO LPV AFIN (aumentado con integrador del acetato)
% estados aumentados:  x = [Z-Z* ; S-S* ; A-A* ; zI],  dzI = -(A - A*)
Cz = [0 0 1];                 % salida: acetato
A0_lin = diag([-alpha*Dop, -Dop, -Dop]);
Ar1 = [ 1 0 0; -1 0 0;  k 0 0];
Ar2 = [ 0 1 0;  0 -1 0; 0 k 0];
Bu  = [-alpha*Zop; Sin-Sop; -k*(Sin-Sop)];   % entrada de control (Ain se cancela)
Bw  = [0; 1; 0];                              % perturbacion (Sin entra en dS)

A0a = [A0_lin, zeros(3,1); -Cz, 0];
A1a = [Ar1, zeros(3,1); zeros(1,4)];
A2a = [Ar2, zeros(3,1); zeros(1,4)];
Bua = [Bu; 0];
Bwa = [Bw; 0];

nx=4; nu=1; nw=1; nz=1;

% >>>>>>>>>>>>  PESOS DE DESEMPENO H-inf  (SINTONIZABLES)  <<<<<<<<<<<<
we = 1.0;                      % peso sobre el error de acetato (salida z)
Czp = we*[0 0 1 0];            % z = we*(A - A*)   (penaliza desviacion de acetato)
Dcu = 0;                       % feedthrough control->z   (0: sin penalizar u aqui)
Dnw = 0;                       % feedthrough perturbacion->z
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%% 5) ESTRUCTURA ROLMIP  (A afin en rho1,rho2 ; B, Bw, Cz constantes)
cellA{1}={[0 0], A0a}; cellA{2}={[1 0], A1a}; cellA{3}={[0 1], A2a};
cellB{1}={[0 0], Bua}; cellB{2}={[1 0], zeros(nx,nu)}; cellB{3}={[0 1], zeros(nx,nu)};
cellBn{1}={[0 0], Bwa};cellBn{2}={[1 0], zeros(nx,nw)};cellBn{3}={[0 1], zeros(nx,nw)};
cellCn{1}={[0 0], Czp};cellCn{2}={[1 0], zeros(nz,nx)};cellCn{3}={[0 1], zeros(nz,nx)};
cellDc{1}={[0 0], Dcu};cellDc{2}={[1 0], 0};           cellDc{3}={[0 1], 0};
cellDn{1}={[0 0], Dnw};cellDn{2}={[1 0], 0};           cellDn{3}={[0 1], 0};

Ac = rolmipvar(cellA ,'Ac',lim);
Bc = rolmipvar(cellB ,'Bc',lim);
Bn = rolmipvar(cellBn,'Bn',lim);
Cn = rolmipvar(cellCn,'Cn',lim);
Dc = rolmipvar(cellDc,'Dc',lim);
Dn = rolmipvar(cellDn,'Dn',lim);

% Variables de decision
%  P robusta (comun),  Z dependiente de parametros -> K(rho) programable
P = rolmipvar(nx,nx,'P','sym', [2 2],[0 0]);
Z = rolmipvar(nu,nx,'Z','full',[2 2],[1 1]);
ro = sdpvar(1,1);     % = gamma^2  (norma H-inf al cuadrado)

%% 6) LMIs:  H-infinito (lema real acotado) + Region D
% >>>>>>>>>>>>  PARAMETROS REGION D  (SINTONIZABLES)  <<<<<<<<<<<<
reg_alpha = 0.2;              % tasa minima de decaimiento
reg_r     = 3;               % radio del disco
reg_theta = deg2rad(40);     % angulo de la cuna (amortiguamiento)
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% --- H-infinito ---
Y11 = Ac*P + P*Ac' + Bc*Z + Z'*Bc';
Y21 = Bn';
Y31 = Cn*P + Dc*Z;
Y22 = -eye(nw);
Y32 = Dn;
Y33 = -ro*eye(nz);
SS = [ Y11 , Y21', Y31';
       Y21 , Y22 , Y32';
       Y31 , Y32 , Y33 ];

% --- Region D ---
ct=cos(reg_theta); st=sin(reg_theta);
Rest1 = Ac*P + P'*Ac' + Bc*Z + Z'*Bc' + 2*reg_alpha*P;
Sym   = Ac*P + P'*Ac' + Bc*Z + Z'*Bc';
Asym  = Ac*P - P'*Ac' + Bc*Z - Z'*Bc';
Rest2 = [ st*Sym,                          ct*Asym;
          ct*(-Ac*P+P'*Ac'-Bc*Z+Z'*Bc'),  st*Sym ];
Rest3 = [ -reg_r*P,        P*Ac'+Z'*Bc';
           Ac*P+Bc*Z,     -reg_r*P ];

LMIs = [ ro>=0, P>=1e-6*eye(nx), SS<=0, Rest1<=0, Rest2<=0, Rest3<=0 ];
optimize(LMIs, ro, sdpsettings('solver','sedumi','verbose',1));
pmin = min(checkset(LMIs));

%% 7) RESULTADO Y GANANCIAS DE VERTICE
fprintf('\n===== RESULTADO H-inf + Region D =====\n');
if pmin >= -1e-7
    fprintf('Factible (pmin = %.3e).  gamma Hinf = %.4f\n', pmin, sqrt(value(ro)));
else
    fprintf('INFACTIBLE (pmin = %.3e). Relajar: bajar we, bajar reg_alpha, subir reg_r/theta\n', pmin);
end

% Orden de vertices V11,V21,V12,V22  (coincide con la interpolacion):
r3dummy = [];  % 2 parametros, no hay tercero
beta  = { [1 0],[1 0];     % V11 (r1min,r2min)
          [0 1],[1 0];     % V21 (r1max,r2min)
          [1 0],[0 1];     % V12 (r1min,r2max)
          [0 1],[0 1] };   % V22 (r1max,r2max)
vrho  = [R1(1) R2(1); R1(2) R2(1); R1(1) R2(2); R1(2) R2(2)];
nomK  = {'K11','K21','K12','K22'};
Kv = nan(4,nx);
fprintf('\n--- Ganancias de vertice ---\n');
for i=1:4
    Pi = eval_rolmip(P, beta(i,:), vrho(i,:));
    Zi = eval_rolmip(Z, beta(i,:), vrho(i,:));
    Kv(i,:) = Zi/Pi;
    fprintf('%s (rho1=%.4f rho2=%.4f) = [% .6g % .6g % .6g % .6g];\n', ...
            nomK{i}, vrho(i,1), vrho(i,2), Kv(i,:));
end

%% 8) REGRESION AFIN:  K(rho) = K0 + K1*rho1 + K2*rho2   (forma para el PLC)
Phi = [ones(4,1), vrho(:,1), vrho(:,2)];     % 4x3
K0=zeros(1,nx); K1=zeros(1,nx); K2=zeros(1,nx);
for j=1:nx
    c = Phi \ Kv(:,j);     % [K0j; K1j; K2j]
    K0(j)=c(1); K1(j)=c(2); K2(j)=c(3);
end
% error de reconstruccion (bilineal vs afin)
err = 0;
for i=1:4
    Krec = K0 + K1*vrho(i,1) + K2*vrho(i,2);
    err = max(err, max(abs(Krec - Kv(i,:))));
end
fprintf('\n--- Gain scheduling afin  K(rho)=K0+K1*rho1+K2*rho2 ---\n');
fprintf('K0 = [% .6g % .6g % .6g % .6g];\n', K0);
fprintf('K1 = [% .6g % .6g % .6g % .6g];\n', K1);
fprintf('K2 = [% .6g % .6g % .6g % .6g];\n', K2);
fprintf('Error de reconstruccion afin = %.2e\n', err);

% guardar para el lazo / el PLC
K11=Kv(1,:); K21=Kv(2,:); K12=Kv(3,:); K22=Kv(4,:);
save('lpv_gains.mat','K11','K21','K12','K22','K0','K1','K2', ...
     'R1','R2','Dop');
fprintf('\n>> Guardado lpv_gains.mat (R1,R2, ganancias de vertice y K0,K1,K2)\n');

% ====================================================================
%  HELPERS (no dependen de evalpar)
% ====================================================================
function M = eval_rolmip(var, beta, rho_phys)
    if exist('evalpar','file')==2
        try
            M = double(evalpar(var, rho_phys));
            if isnumeric(M) && ~any(isnan(M(:))), return; end
        catch
        end
    end
    D = rolmip_data(var);
    M = 0;
    for m = 1:numel(D)
        e = D(m).exponent;
        if ~iscell(e)
            tmp = cell(1,numel(beta)); p = 1;
            for s = 1:numel(beta), g = numel(beta{s}); tmp{s}=e(p:p+g-1); p=p+g; end
            e = tmp;
        end
        w = 1;
        for s = 1:numel(beta), w = w*prod(beta{s}.^e{s}); end
        if w==0, continue; end
        c = D(m).value; if isa(c,'sdpvar'), c = value(c); end
        M = M + w*c;
    end
end
function D = rolmip_data(var)
    try, D = var.data; return; catch, end
    try, s = struct(var); D = s.data; return; catch, end
    error('No pude leer var.data en tu ROLMIP. Ejecuta fieldnames(P) y pasamelo.');
end
