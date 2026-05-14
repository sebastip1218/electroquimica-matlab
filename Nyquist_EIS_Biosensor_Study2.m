%% Nyquist_EIS_Biosensor_Study.m
% Estudio visual de EIS para un biosensor tipo Randles:
% Rs + (Rct || Cdl) + Warburg
%
% Qué muestra:
% 1) Nyquist SIN Warburg (para ver limpio Rs y Rct)
% 2) Nyquist CON Warburg (para ver la cola difusiva)
% 3) Bode magnitud
% 4) Bode fase
% 5) Curva de calibración DeltaRct vs log C
%
% IMPORTANTE:
% - En Nyquist:
%     eje x = Z' (ohm)
%     eje y = -Z'' (ohm)
% - En Bode:
%     eje x = frecuencia (Hz)
%     eje y = magnitud o fase (grados)
%
% Cambia SOLO el bloque de parámetros de abajo para explorar el sistema.

clear; clc; close all;

%% =========================
%  BLOQUE DE PARÁMETROS
%  (MODIFICA AQUÍ)
%  =========================

% Electrode and electrochemical system
r_elec = 1.5e-3;       % m, radio del electrodo (3 mm de diametro => 1.5 mm de radio)
A = pi*r_elec^2;       % m^2, area geometrica
Rs = 60;               % ohm, resistencia de solucion
Cdl0 = 30e-6;          % F, capacitancia de doble capa base
sigmaW = 300;          % ohm*s^(-1/2), parametro Warburg
Vac = 5e-3;            % V, perturbacion AC (5 mV)
eta_dc = 0.0;          % V, sobrepotencial DC alrededor del cual se perturba

% Frequency window for EIS
fmin = 1e-4;           % Hz
fmax = 1e4;            % Hz
Nf = 400;              % numero de puntos de frecuencia

% Kinetic + redox probe parameters
k0 = 2e-5;             % m/s, velocidad estandar de transferencia de carga
Credox_mM = 5;         % mM, concentracion del par ferri/ferro

% Butler-Volmer (solo como valor ilustrativo)
alpha = 0.5;           % coeficiente de transferencia
eta_demo = 50e-3;      % V, sobrepotencial de ejemplo para calcular i_BV

% Biosensor blocking model (modelo educativo)
Kd_target = 5e-13;     % M, "afinidad" aparente del target por la sonda
beta_block = 25;       % cuanto crece Rct por bloqueo / hibridacion
gamma_Cdl = 0.30;      % fraccion con la que disminuye Cdl al bloquearse

% Concentraciones del biomarcador objetivo (M)
Ctarget_list = [0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11];
Ctarget_labels = {'Blank', '1 fM', '10 fM', '100 fM', '1 pM', '10 pM'};

%% =========================
%  CONSTANTES
%  =========================
R = 8.314462618;       % J mol^-1 K^-1
T = 298.15;            % K
F = 96485.33212;       % C mol^-1
n = 1;                 % numero de electrones

%% =========================
%  FRECUENCIAS
%  =========================
% Se usa orden descendente: alta frecuencia -> baja frecuencia
f = logspace(log10(fmax), log10(fmin), Nf);
omega = 2*pi*f;

%% =========================
%  CALCULO DE i0 Y Rct BASE
%  =========================
% 1 mM = 1 mol/m^3
Credox = Credox_mM;  % mol/m^3

% Corriente de intercambio educativa:
% i0 = n F A k0 C*
i0 = n*F*A*k0*Credox;              % A

% Resistencia de transferencia de carga base:
% Rct = RT / (n F i0)
Rct0 = R*T/(n*F*i0);               % ohm

% Corriente Butler-Volmer de ejemplo
i_BV_demo = i0*( exp((1-alpha)*n*F*eta_demo/(R*T)) - exp(-alpha*n*F*eta_demo/(R*T)) );

%% =========================
%  PREALOCACION
%  =========================
nCases = numel(Ctarget_list);
Z_noW_all = cell(nCases,1);
Z_withW_all = cell(nCases,1);

results = table('Size',[nCases, 10], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Label','Ctarget_M','Theta','Rct0_ohm','RctEff_ohm','CdlEff_uF','Tau_s','fPeak_Hz','i0_uA','iBV_demo_uA'});

%% =========================
%  BUCLE PRINCIPAL
%  =========================
for k = 1:nCases
    
    Ct = Ctarget_list(k);
    
    % Fraccion de cobertura / bloqueo (modelo educativo)
    if Ct == 0
        theta = 0;
    else
        theta = Ct/(Ct + Kd_target);
    end
    
    % Rct efectiva: aumenta al aumentar el bloqueo
    Rct_eff = Rct0*(1 + beta_block*theta);
    
    % Cdl efectiva: disminuye al aumentar el bloqueo
    Cdl_eff = Cdl0*(1 - gamma_Cdl*theta);
    Cdl_eff = max(Cdl_eff, 1e-9);   % evita valores no fisicos
    
    % Constante de tiempo y frecuencia aproximada del maximo imaginario
    tau = Rct_eff*Cdl_eff;
    fPeak = 1/(2*pi*tau);
    
    % Impedancia capacitiva ideal
    Zc = 1./(1i*omega*Cdl_eff);
    
    % Paralelo Rct || Cdl
    Zpar = 1./(1./Rct_eff + 1./Zc);
    
    % Warburg semi-infinito
    Zw = sigmaW*(1 - 1i)./sqrt(omega);
    
    % Impedancia total
    Z_noW = Rs + Zpar;
    Z_withW = Rs + Zpar + Zw;
    
    % Guardar curvas
    Z_noW_all{k} = Z_noW;
    Z_withW_all{k} = Z_withW;
    
    % Guardar resultados en tabla (CORREGIDO)
    results.Label(k)      = string(Ctarget_labels{k});
    results.Ctarget_M(k)  = Ct;
    results.Theta(k)      = theta;
    results.Rct0_ohm(k)   = Rct0;
    results.RctEff_ohm(k) = Rct_eff;
    results.CdlEff_uF(k)  = Cdl_eff*1e6;
    results.Tau_s(k)      = tau;
    results.fPeak_Hz(k)   = fPeak;
    results.i0_uA(k)      = i0*1e6;
    results.iBV_demo_uA(k)= i_BV_demo*1e6;
end

%% =========================
%  ELEGIR UNA CURVA PARA ANOTAR
%  =========================
idxAnnot = nCases;     % usa la ultima concentracion para anotar
Zann_noW = Z_noW_all{idxAnnot};
Zann_W   = Z_withW_all{idxAnnot};

% Punto de maximo imaginario
[~, idxMaxImag] = max(-imag(Zann_noW));
Zx_peak = real(Zann_noW(idxMaxImag));
Zy_peak = -imag(Zann_noW(idxMaxImag));

%% =========================
%  GRAFICAS PRINCIPALES
%  =========================
colors = lines(nCases);

fig1 = figure('Color','w','Name','EIS Study - Nyquist and Bode');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% ---------------------------------
% 1) Nyquist SIN Warburg
% ---------------------------------
nexttile;
hold on; box on; grid on;
for k = 1:nCases
    Zp = Z_noW_all{k};
    plot(real(Zp), -imag(Zp), 'LineWidth', 1.8, 'Color', colors(k,:));
end
xlabel('Z'' (\Omega)');
ylabel('-Z'''' (\Omega)');
title('Nyquist SIN Warburg (lectura limpia de R_s y R_{ct})');

% Anotaciones
plot(real(Zann_noW(1)), -imag(Zann_noW(1)), 'ko', 'MarkerFaceColor','k');
text(real(Zann_noW(1))*1.02, 0.03*max(-imag(Zann_noW)), ...
    'Alta f \rightarrow R_s', 'FontSize', 10, 'Color','k');

plot(Zx_peak, Zy_peak, 'ko', 'MarkerFaceColor','y');
text(Zx_peak*1.02, Zy_peak*1.02, ...
    'Punto de -Z'''' máximo', 'FontSize', 10, 'Color','k');

% Anotar Rct horizontal
x_left  = real(Zann_noW(1));
x_right = real(Zann_noW(end));
y_arrow = -0.08*max(-imag(Zann_noW));
plot([x_left x_right], [0 0], 'k.', 'MarkerSize', 14);
annotation('doublearrow', [0.15 0.32], [0.49 0.49]); %#ok<UNRCH> 
% La annotation depende del layout de la figura; por eso también ponemos texto dentro del eje:
text((x_left+x_right)/2, 0.07*max(-imag(Zann_noW)), ...
    sprintf('Diametro horizontal \\approx R_{ct} = %.1f \\Omega', results.RctEff_ohm(idxAnnot)), ...
    'HorizontalAlignment','center', 'FontSize', 10, 'BackgroundColor','w');

legend(Ctarget_labels, 'Location','best');

% ---------------------------------
% 2) Nyquist CON Warburg
% ---------------------------------
nexttile;
hold on; box on; grid on;
for k = 1:nCases
    Zp = Z_withW_all{k};
    plot(real(Zp), -imag(Zp), 'LineWidth', 1.8, 'Color', colors(k,:));
end
xlabel('Z'' (\Omega)');
ylabel('-Z'''' (\Omega)');
title('Nyquist CON Warburg (cola difusiva a baja frecuencia)');

plot(real(Zann_W(1)), -imag(Zann_W(1)), 'ko', 'MarkerFaceColor','k');
text(real(Zann_W(1))*1.02, 0.03*max(-imag(Zann_W)), ...
    'Alta f', 'FontSize', 10, 'Color','k');

plot(real(Zann_W(end)), -imag(Zann_W(end)), 'ko', 'MarkerFaceColor','c');
text(real(Zann_W(end))*0.92, -imag(Zann_W(end))*0.92, ...
    'Baja f \rightarrow difusion / Warburg', 'FontSize', 10, 'Color','k');

legend(Ctarget_labels, 'Location','best');

% ---------------------------------
% 3) Bode magnitud (para una curva)
% ---------------------------------
nexttile;
hold on; box on; grid on;
semilogx(f, abs(Zann_W), 'b-', 'LineWidth', 2);
xlabel('Frecuencia (Hz)');
ylabel('|Z| (\Omega)');
title(['Bode magnitud - ' char(results.Label(idxAnnot))]);

% Marcar frecuencia caracteristica aproximada
xline(results.fPeak_Hz(idxAnnot), '--r', 'f_{peak}', ...
    'LabelOrientation','horizontal', 'LabelVerticalAlignment','middle');

% ---------------------------------
% 4) Bode fase (para una curva)
% ---------------------------------
nexttile;
hold on; box on; grid on;
phase_deg = angle(Zann_W)*180/pi;
semilogx(f, phase_deg, 'm-', 'LineWidth', 2);
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
title(['Bode fase - ' char(results.Label(idxAnnot))]);

xline(results.fPeak_Hz(idxAnnot), '--r', 'f_{peak}', ...
    'LabelOrientation','horizontal', 'LabelVerticalAlignment','middle');

%% =========================
%  CURVA DE CALIBRACION
%  =========================
% Señal analitica educativa:
% DeltaRct = Rct( muestra ) - Rct( blanco )
Rct_blank = results.RctEff_ohm(1);
DeltaRct = results.RctEff_ohm - Rct_blank;

% Excluir el blanco para el ajuste logaritmico
valid = results.Ctarget_M > 0;
C_pM = results.Ctarget_M(valid)*1e12;    % convertir M -> pM
logC_pM = log10(C_pM);
DeltaRct_kohm = DeltaRct(valid)/1000;

% Ajuste lineal: DeltaRct = m*logC + b
p = polyfit(logC_pM, DeltaRct_kohm, 1);
m = p(1);
b = p(2);

xfit = linspace(min(logC_pM), max(logC_pM), 200);
yfit = polyval(p, xfit);

% Calculo de R^2
ycalc = polyval(p, logC_pM);
SSres = sum((DeltaRct_kohm - ycalc).^2);
SStot = sum((DeltaRct_kohm - mean(DeltaRct_kohm)).^2);
R2 = 1 - SSres/SStot;

fig2 = figure('Color','w','Name','Calibration Curve');
hold on; box on; grid on;
plot(logC_pM, DeltaRct_kohm, 'ko', 'MarkerFaceColor','g', 'MarkerSize', 8);
plot(xfit, yfit, 'r-', 'LineWidth', 2);
xlabel('log C_{target} (pM)');
ylabel('\DeltaR_{ct} (k\Omega)');
title('\DeltaR_{ct} vs log C (modelo educativo)');

for k = find(valid)'
    text(log10(results.Ctarget_M(k)*1e12), DeltaRct(k)/1000, ['  ' char(results.Label(k))], 'FontSize', 9);
end

txtEq = sprintf('\\DeltaR_{ct} = %.3f logC + %.3f\nR^2 = %.4f', m, b, R2);
text(mean(xfit), max(yfit)*0.85, txtEq, 'FontSize', 11, 'BackgroundColor','w');

%% =========================
%  MOSTRAR TABLA DE RESULTADOS
%  =========================
disp(' ');
disp('================ RESULTADOS PRINCIPALES ================');
disp(results);

disp(' ');
disp('================ FORMULAS USADAS ================');
disp('i0 = n F A k0 C*');
disp('Rct0 = RT / (n F i0)');
disp('i_BV = i0 [exp((1-alpha)nFeta/RT) - exp(-alpha nFeta/RT)]');
disp('Zc = 1 / (j w Cdl)');
disp('Zw = sigma (1 - j) / sqrt(w)');
disp('Ztotal = Rs + (Rct || Cdl) + Zw');
disp(' ');
disp('================ INTERPRETACION RAPIDA ================');
disp('- Mayor k0  -> menor Rct -> semicírculo más pequeño');
disp('- Mayor concentracion del biomarcador -> mayor bloqueo -> mayor Rct -> semicírculo más grande');
disp('- Mayor sigmaW o menor fmin -> cola de Warburg más visible');
disp('- En Nyquist: izquierda = alta frecuencia; derecha = baja frecuencia');
disp('- En Bode: eje x = frecuencia (Hz); eje y = magnitud o fase');

%% =========================
%  TABLA EXTRA RESUMIDA
%  =========================
summaryTable = table(results.Label, results.Ctarget_M*1e12, results.Theta, ...
    results.RctEff_ohm, DeltaRct, results.CdlEff_uF, results.fPeak_Hz, ...
    'VariableNames', {'Label','Ctarget_pM','Theta','RctEff_ohm','DeltaRct_ohm','CdlEff_uF','fPeak_Hz'});

disp(' ');
disp('================ TABLA RESUMIDA ================');
disp(summaryTable);
disp(' ');
fprintf('Pendiente m = %.4f kOhm/decada\n', m);
fprintf('Intercepto b = %.4f kOhm\n', b);
fprintf('R^2 = %.4f\n', R2);

%% =========================
%  NOTAS IMPORTANTES
%  =========================
% 1) En Nyquist SIN Warburg:
%    - Intercepto izquierdo ~ Rs
%    - Diametro horizontal del semicírculo ~ Rct
%
% 2) En Nyquist CON Warburg:
%    - A bajas frecuencias aparece la cola difusiva
%    - Ya no siempre se "lee" tan limpio Rct visualmente
%
% 3) En un biosensor real:
%    - El biomarcador no entra directamente en la formula de Butler-Volmer
%    - El biomarcador bloquea la interfase y por eso cambia Rct
%
% 4) Este script es educativo.
%    - Ayuda a entender tendencias
%    - No sustituye un ajuste riguroso experimental
