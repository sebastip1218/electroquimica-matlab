%% Tarea 2 - Polarización mixta (BV + transporte de masa) - Código maestro
% Autor: (tu nombre)
% Objetivo:
% - Calcula i0 desde j0*A
% - Convierte concentraciones M -> mol/cm^3
% - Calcula corrientes límite i_l,c e i_l,a con i_l = n*F*A*m*C*
% - Evalúa i(E) con BV + transporte de masa (convención IUPAC)
% - Grafica cada contribución y la corriente total

clear; clc; close all;

%% ------------------------- 1) Constantes -------------------------
F = 96485;            % C/mol
R = 8.314;            % J/(mol*K)
T = 298.15;           % K

alpha = 0.5;          % dado en la tarea
A = 0.1;              % cm^2 (área del electrodo)
mO = 0.01;            % cm/s (transporte de masa para especie Ox)
mR = 0.01;            % cm/s (transporte de masa para especie Red)

tinyI = 1e-12;        % A  (para "corriente límite ~0" sin dividir entre 0)

%% ------------------------- 2) Dominio de potencial E -------------------------
% Ajusta el rango y el paso según tu tarea (puntos from-to)
Emin = -2.0;          % V
Emax = +2.0;          % V
dE   =  0.01;         % V
E = (Emin:dE:Emax).'; % columna
N = numel(E);
fprintf('Puntos en el eje X: %d\n', N);

%% ------------------------- 3) Composición de la solución (M) -------------------------
% Solución: 1.0 M HBr + 0.1 M K3Fe(CN)6
% (asumiendo disociación completa para HBr)
C_Hplus_M   = 1.0;    % M
C_H2_M      = 0.0;    % M (no hay H2 al inicio)
C_Brminus_M = 1.0;    % M
C_Br2_M     = 0.0;    % M (no hay Br2 al inicio)
C_Fe3_M     = 0.1;    % M  Fe(CN)6^3- (Ox)
C_Fe4_M     = 0.0;    % M  Fe(CN)6^4- (Red)

% Conversión M -> mol/cm^3
M_to_molcm3 = 1/1000; % 1 M = 1e-3 mol/cm^3
C_Hplus   = C_Hplus_M   * M_to_molcm3;
C_H2      = C_H2_M      * M_to_molcm3;
C_Brminus = C_Brminus_M * M_to_molcm3;
C_Br2     = C_Br2_M     * M_to_molcm3;
C_Fe3     = C_Fe3_M     * M_to_molcm3;
C_Fe4     = C_Fe4_M     * M_to_molcm3;

%% ------------------------- 4) Definir parejas redox (parámetros) -------------------------
% IMPORTANTE: Para que tus curvas se desplacen como en tu figura,
% en tus chats usaron "potenciales de oxidación" (negativos):
%   H:  Eeq =  0.000 V
%   Br: Eeq = + 1.087 V
%   Fe: Eeq = + 0.361 V
%
% Si tu profe quiere potenciales de REDUCCIÓN (positivos), cámbialos a:
%   Br: +1.087, Fe: +0.361 (H sigue en 0)
%
% Además:
% - nKinetics: n usado en f = nF/RT dentro de BV
% - nLim: n usado en i_l = nF A m C* (en tu tarea lo usaste consistente con cada pareja)
%
% j0 se da en A/cm^2 (densidad), luego i0 = j0*A

rxn(1).name = 'H^+/H_2';
rxn(1).nKinetics = 2;                  % n para f (como usaste en tus gráficas)
rxn(1).nLim      = 2;                  % n para i_l (si quieres 0.19297 A con 1 M)
rxn(1).j0 = 1e-3;                      % A/cm^2  (dado)
rxn(1).Eeq = 0.0;                      % V (oxidación = 0)
rxn(1).C_Ox = C_Hplus;                 % Ox disponible para reducción (H+)
rxn(1).C_Red = C_H2;                   % Red disponible para oxidación (H2 ~ 0)

rxn(2).name = 'Br_2/Br^-';
rxn(2).nKinetics = 2;
rxn(2).nLim      = 2;
rxn(2).j0 = 1e-2;                      % A/cm^2 (dado)
rxn(2).Eeq = +1.087;                   % V (oxidación, como tu figura)
rxn(2).C_Ox = C_Br2;                   % Br2 ~ 0
rxn(2).C_Red = C_Brminus;              % Br- = 1 M

rxn(3).name = 'Fe(CN)_6^{3-}/Fe(CN)_6^{4-}';
rxn(3).nKinetics = 1;
rxn(3).nLim      = 1;
rxn(3).j0 = 4e-5;                      % A/cm^2 (dado)
rxn(3).Eeq = +0.361;                   % V (oxidación, como tu figura)
rxn(3).C_Ox = C_Fe3;                   % 0.1 M
rxn(3).C_Red = C_Fe4;                  % ~0 M

% (Opcional) si quieres usar Nernst cuando existan ambas especies:
useNernst = false;  % pon true si te lo piden y tienes ambas especies no-cero

%% ------------------------- 5) Cálculo de i(E) para cada reacción -------------------------
I = zeros(N, numel(rxn));  % cada columna será una reacción

disp([rxn.Eeq])
for k = 1:numel(rxn)
    % 5.1 i0 desde j0*A
    i0 = rxn(k).j0 * A;

    % 5.2 corrientes límite desde i_l = n F A m C*
    % Convención IUPAC:
    % - i_l,c (reducción) NEGATIVA, depende de Ox (reactivo en reducción)
    % - i_l,a (oxidación) POSITIVA, depende de Red (reactivo en oxidación)
    ilc_mag = rxn(k).nLim * F * A * mO * rxn(k).C_Ox;  % A (magnitud)
    ila_mag = rxn(k).nLim * F * A * mR * rxn(k).C_Red; % A (magnitud)

    ilc = -max(ilc_mag, tinyI);  % negativo
    ila = +max(ila_mag, tinyI);  % positivo

    % 5.3 Eeq: o fijo (de la tabla) o Nernst si procede
    Eeq = rxn(k).Eeq;
    if useNernst
        % Para Ox + ne- -> Red:
        % Eeq = E0 - (RT/nF)*ln(a_Red/a_Ox)
        if rxn(k).C_Ox > 0 && rxn(k).C_Red > 0
            Q = rxn(k).C_Red / rxn(k).C_Ox;
            Eeq = rxn(k).Eeq - (R*T/(rxn(k).nKinetics*F))*log(Q);
        end
    end

    % 5.4 sobrepotencial eta = E - Eeq
    eta = E - Eeq;

    % 5.5 f = nF/RT
    f = rxn(k).nKinetics * F / (R*T);

    % 5.6 i(E) con BV + transporte de masa (forma estable)
    I(:,k) = bv_mass_transport_current(eta, i0, alpha, f, ila, ilc);
    
    % imprimir resumen
    fprintf('\n%s\n', rxn(k).name);
    fprintf('  i0  = %.4g A\n', i0);
    fprintf('  ila = %.4g A\n', ila);
    fprintf('  ilc = %.4g A\n', ilc);
    fprintf('  Eeq = %.4g V\n', Eeq);
end

%% ------------------------- 6) Corriente total -------------------------
I_total = sum(I, 2);

%% ------------------------- 7) Gráficas -------------------------
figure('Color','w');
hold on; grid on;

% Graficar individuales
for k = 1:numel(rxn)
    plot(E, I(:,k), 'LineWidth', 1.6);
end

% Total (línea punteada)
plot(E, I_total, '--', 'LineWidth', 2.2);

xlabel('Potencial del electrodo, E (V)');
ylabel('Corriente, i (A)');
title('Polarización mixta: contribuciones individuales y corriente total');

leg = cell(1, numel(rxn)+1);
for k = 1:numel(rxn), leg{k} = rxn(k).name; end
leg{end} = 'Total';
legend(leg, 'Location','best');

%% (Opcional) Guardar figura
% saveas(gcf, 'tarea2_polarizacion_mixta.png');

%% ------------------------- Función local -------------------------
function i = bv_mass_transport_current(eta, i0, alpha, f, ila, ilc)
% BV con transporte de masa:
% i = [exp(-a f eta) - exp((1-a) f eta)] / [1/i0 + exp(-a f eta)/ilc - exp((1-a) f eta)/ila]
%
% Convención IUPAC requerida:
% - ila > 0
% - ilc < 0

    % Evitar overflow en exp() (seguro numérico)
    x1 = -alpha*f.*eta;
    x2 = (1-alpha)*f.*eta;
    x1 = clamp(x1, -700, 700);
    x2 = clamp(x2, -700, 700);

    num = exp(x1) - exp(x2);
    den = (1./i0) + (exp(x1)./ilc) - (exp(x2)./ila);

    i = num ./ den;
end

function y = clamp(x, lo, hi)
    y = min(max(x, lo), hi);
end
