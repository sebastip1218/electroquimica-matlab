%% DOE_Hu_MOF_EIS_Simulator.m
% Simulador educativo de EIS para un biosensor CoNi-MOF/cDNA/miRNA
% Basado en el esquema de Hu et al.: Au electrode + CoNi-MOF(1:1) + cDNA + miRNA.
% Objetivo: generar curvas Nyquist/Bode y una tabla DOE para biomarcador, concentración y ventana de frecuencia.

clear; clc; close all;
rng(8);

%% =========================================================
%  1) CONDICIONES FIJAS DEL SENSOR
%  =========================================================
% Estas condiciones se mantienen constantes en el DOE principal.
% Si se cambian, ya serían factores adicionales.

Rs = 61.0;                  % ohm, resistencia de solución aproximada
Rct_cDNA = 494.5;            % ohm, estado base cDNA/CoNi-MOF/AE antes del target
Cdl0 = 25e-6;                % F, capacitancia efectiva inicial
sigmaW = 250;                % ohm*s^(-1/2), Warburg educativo
Nf = 400;                    % puntos por barrido
nRep = 3;                    % réplicas por corrida
noise_Rct_CV = 0.015;        % ruido relativo experimental simulado en Rct
noise_Cdl_CV = 0.010;        % ruido relativo experimental simulado en Cdl

%% =========================================================
%  2) FACTORES DEL DOE
%  =========================================================
% Factor A: tipo de biomarcador.
% El cDNA se considera diseñado para miRNA-126; los otros biomarcadores se simulan
% como interferentes/similares con menor afinidad y menor bloqueo.

biomarker = table;
biomarker.Name = ["miRNA-126"; "miRNA-21"; "miRNA-155"];
biomarker.Kd_fM = [1000; 5000; 8000];          % afinidad aparente: menor Kd = mayor unión
biomarker.DeltaRctMax_kOhm = [5.60; 1.00; 0.70]; % bloqueo máximo posible
biomarker.CdlDropMax = [0.30; 0.12; 0.08];     % caída máxima relativa de Cdl

% Factor B: concentración.
% Se eligen puntos en escala lineal y en una ventana estrecha para aproximar respuesta lineal.
C_fM_levels = [10 30 50 70 90];

% Factor C: protocolo de frecuencia.
% Hu completo replica el rango del paper. El recortado sirve para probar si se pierde información.
freqWindow = table;
freqWindow.Name = ["Hu completo"; "Sin ultra-baja f"];
freqWindow.fmin_Hz = [1e-2; 1];
freqWindow.fmax_Hz = [1e5; 1e5];

%% =========================================================
%  3) CONSTRUIR MATRIZ DOE
%  =========================================================
row = 0;
DOE = table;
for a = 1:height(biomarker)
    for b = 1:numel(C_fM_levels)
        for c = 1:height(freqWindow)
            for r = 1:nRep
                row = row + 1;
                DOE.Run(row,1) = row;
                DOE.Biomarker(row,1) = biomarker.Name(a);
                DOE.Concentration_fM(row,1) = C_fM_levels(b);
                DOE.FrequencyProtocol(row,1) = freqWindow.Name(c);
                DOE.fmin_Hz(row,1) = freqWindow.fmin_Hz(c);
                DOE.fmax_Hz(row,1) = freqWindow.fmax_Hz(c);
                DOE.Replica(row,1) = r;
            end
        end
    end
end

%% =========================================================
%  4) SIMULAR EIS PARA CADA CORRIDA
%  =========================================================
Zcurves = cell(height(DOE),1);
FreqCurves = cell(height(DOE),1);

for i = 1:height(DOE)
    a = find(biomarker.Name == DOE.Biomarker(i));
    Ct_fM = DOE.Concentration_fM(i);

    % Modelo de unión tipo Langmuir, usado en una región donde se aproxima a lineal.
    theta = Ct_fM/(Ct_fM + biomarker.Kd_fM(a));

    DeltaRct_true = 1000 * biomarker.DeltaRctMax_kOhm(a) * theta;
    Rct_true = Rct_cDNA + DeltaRct_true;
    Cdl_true = Cdl0 * (1 - biomarker.CdlDropMax(a)*theta);

    % Ruido por réplica, para que las réplicas no sean idénticas.
    Rct_obs = Rct_true * (1 + noise_Rct_CV*randn);
    Cdl_obs = Cdl_true * (1 + noise_Cdl_CV*randn);
    Cdl_obs = max(Cdl_obs, 1e-9);

    f = logspace(log10(DOE.fmax_Hz(i)), log10(DOE.fmin_Hz(i)), Nf);
    w = 2*pi*f;

    Zc = 1./(1i*w*Cdl_obs);
    Zpar = 1./(1./Rct_obs + 1./Zc);
    Zw = sigmaW*(1 - 1i)./sqrt(w);
    Z = Rs + Zpar + Zw;

    % Lectura visual aproximada del diámetro observado en la ventana.
    % Ojo: no sustituye el ajuste de circuito equivalente.
    Rct_visual = max(real(Z)) - min(real(Z));
    DeltaRct_visual = Rct_visual - Rct_cDNA;

    tau = Rct_obs*Cdl_obs;
    fPeak = 1/(2*pi*tau);

    DOE.Theta(i,1) = theta;
    DOE.Rct_true_ohm(i,1) = Rct_true;
    DOE.Rct_obs_ohm(i,1) = Rct_obs;
    DOE.DeltaRct_obs_ohm(i,1) = Rct_obs - Rct_cDNA;
    DOE.Cdl_obs_uF(i,1) = Cdl_obs*1e6;
    DOE.fPeak_Hz(i,1) = fPeak;
    DOE.Rct_visual_ohm(i,1) = Rct_visual;
    DOE.DeltaRct_visual_ohm(i,1) = DeltaRct_visual;

    Zcurves{i} = Z;
    FreqCurves{i} = f;
end

%% =========================================================
%  5) TABLAS RESUMEN PARA DOE/ANOVA
%  =========================================================
summaryDOE = groupsummary(DOE, {'Biomarker','Concentration_fM','FrequencyProtocol'}, ...
    {'mean','std'}, {'DeltaRct_obs_ohm','Rct_obs_ohm','Cdl_obs_uF','fPeak_Hz','DeltaRct_visual_ohm'});

disp('================ MATRIZ DOE ================');
disp(DOE(:, {'Run','Biomarker','Concentration_fM','FrequencyProtocol','Replica','DeltaRct_obs_ohm','Rct_obs_ohm','Cdl_obs_uF','fPeak_Hz'}));

disp('================ RESUMEN POR CORRIDA ================');
disp(summaryDOE);

writetable(DOE, 'DOE_Hu_MOF_EIS_raw_runs.csv');
writetable(summaryDOE, 'DOE_Hu_MOF_EIS_summary.csv');

%% =========================================================
%  6) GRÁFICA 1: NYQUIST PARA COMPARAR BIOMARCADORES
%  =========================================================
figure('Color','w','Name','Nyquist biomarcadores');
hold on; box on; grid on;
plotRows = find(DOE.Concentration_fM == 50 & DOE.FrequencyProtocol == "Hu completo" & DOE.Replica == 1);
colors = lines(numel(plotRows));
for k = 1:numel(plotRows)
    idx = plotRows(k);
    Z = Zcurves{idx};
    plot(real(Z)/1000, -imag(Z)/1000, 'LineWidth', 2, 'Color', colors(k,:));
end
xlabel('Z'' (k\Omega)');
ylabel('-Z'''' (k\Omega)');
title('EIS simulado: efecto del tipo de biomarcador a 50 fM');
legend(cellstr(DOE.Biomarker(plotRows)), 'Location','best');

%% =========================================================
%  7) GRÁFICA 2: CURVA DE CALIBRACIÓN LINEAL APROXIMADA
%  =========================================================
figure('Color','w','Name','Calibracion DOE');
hold on; box on; grid on;
for a = 1:height(biomarker)
    tmp = summaryDOE(summaryDOE.Biomarker == biomarker.Name(a) & summaryDOE.FrequencyProtocol == "Hu completo", :);
    x = tmp.Concentration_fM;
    y = tmp.mean_DeltaRct_obs_ohm/1000;
    plot(x, y, '-o', 'LineWidth', 2, 'MarkerFaceColor','w');

    p = polyfit(x, y, 1);
    yfit = polyval(p, x);
    SSres = sum((y - yfit).^2);
    SStot = sum((y - mean(y)).^2);
    R2 = 1 - SSres/SStot;
    fprintf('%s: DeltaRct(kOhm) = %.5f*C(fM) + %.5f, R2 = %.4f\n', biomarker.Name(a), p(1), p(2), R2);
end
xlabel('Concentración del biomarcador (fM)');
ylabel('\DeltaR_{ct} observado (k\Omega)');
title('Calibración simulada en escala lineal');
legend(cellstr(biomarker.Name), 'Location','northwest');

%% =========================================================
%  8) GRÁFICA 3: EFECTO DE LA VENTANA DE FRECUENCIA
%  =========================================================
figure('Color','w','Name','Frecuencia como protocolo');
hold on; box on; grid on;
rowsFull = find(DOE.Biomarker == "miRNA-126" & DOE.Concentration_fM == 70 & DOE.Replica == 1);
for k = 1:numel(rowsFull)
    idx = rowsFull(k);
    Z = Zcurves{idx};
    plot(real(Z)/1000, -imag(Z)/1000, 'LineWidth', 2);
end
xlabel('Z'' (k\Omega)');
ylabel('-Z'''' (k\Omega)');
title('Efecto de cambiar el barrido de frecuencias');
legend(cellstr(DOE.FrequencyProtocol(rowsFull)), 'Location','best');

%% =========================================================
%  9) NOTA PARA MINITAB
%  =========================================================
disp(' ');
disp('Para Minitab/ANOVA: usar como Y principal DeltaRct_obs_ohm.');
disp('Factores: Biomarker, Concentration_fM, FrequencyProtocol.');
disp('Si FrequencyProtocol no se quiere como factor principal, mantener Hu completo y correr ANOVA de dos factores.');
