%% Nyquist / Bode educational model for a CoNi-MOF biosensor
% This script is EDUCATIONAL. It is meant to help visualize:
%   1) Rs at high frequency
%   2) Rct + Cdl/CPE semicircle behavior
%   3) Warburg diffusion tail at low frequency
%   4) How k0, redox-probe concentration and target concentration affect the Nyquist plot
%
% IMPORTANT CONCEPTS
% Nyquist axes:
%   x-axis  -> Z'      (real part of impedance, resistive contribution)
%   y-axis  -> -Z''    (negative imaginary part, capacitive/reactive contribution)
%
% Bode axes:
%   x-axis  -> frequency (Hz)
%   y-axis  -> |Z| (ohm) or phase (degrees)
%
% Main electrochemical relations used:
%   Butler-Volmer current (general):
%      i = i0*[exp((1-alpha)*n*F*eta/(R*T)) - exp(-alpha*n*F*eta/(R*T))]
%
%   Exchange current (simple educational estimate for a reversible redox probe):
%      i0 = n*F*A*k0*C*
%
%   Charge-transfer resistance near equilibrium:
%      Rct = R*T / (n*F*i0) = R*T / (n^2*F^2*A*k0*C*)
%
%   Capacitive current:
%      iC = C * dV/dt
%
%   Capacitor impedance:
%      Zc = 1/(j*w*C)
%
%   Warburg impedance (semi-infinite diffusion):
%      Zw = sigma*(1 - j)/sqrt(w)
%
%   Time constant of the RC interfacial process:
%      tau = Rct*Cdl
%      f_peak ~ 1/(2*pi*Rct*Cdl)   (ideal semicircle, no Warburg)
%
% NOTE:
%   The biosensor calibration in the Hu paper is based on the increase of Rct after
%   hybridization, not on ferri/ferro concentration itself. To mimic that idea,
%   this script uses a simple blocking model:
%      theta = Ctarget/(Ctarget + Kd)
%      Rct_eff = Rct0*(1 + beta_block*theta)
%      Cdl_eff = Cdl0*(1 - gamma_Cdl*theta)
%
%   This is not a fitted model to the paper. It is a clean teaching model.

clear; clc; close all;

%% ========================== USER PARAMETERS ==========================
% Universal constants
R  = 8.314;            % J mol^-1 K^-1
F  = 96485;            % C mol^-1
T  = 298.15;           % K
n  = 1;                % number of electrons
alpha = 0.50;          % transfer coefficient

% Electrode and electrochemical system
r_elec = 1.5e-3;       % m, 3 mm diameter Au electrode -> radius = 1.5 mm
A = pi*r_elec^2;       % m^2
Rs = 60;               % ohm, solution resistance (adjustable)
Cdl0 = 12e-6;          % F, double-layer capacitance before target binding
sigmaW = 120;           % ohm*s^( -1/2 )  Warburg coefficient (adjustable)
Vac = 5e-3;            % V, small AC perturbation amplitude used in EIS
eta_dc = 0.0;          % V, DC overpotential for Butler-Volmer current display only

% Frequency window for EIS
fmin = 1e-3;           % Hz
fmax = 1e5;            % Hz
Nf = 400;              % number of frequency points
f = logspace(log10(fmax), log10(fmin), Nf);  % high -> low frequency
w = 2*pi*f;

% Kinetic + redox probe parameters
k0 = 2e-5;             % m s^-1, standard heterogeneous rate constant
Credox_mM = 5;         % mM ferri/ferrocyanide concentration (adjustable)
Credox = Credox_mM;    % keep name for display
Credox_SI = Credox_mM; % mM for display
Credox_M = Credox_mM*1e-3;     % mol L^-1
Credox_m3 = Credox_M*1000;      % mol m^-3

% Biosensor blocking model (educational)
Kd_target = 5e-13;     % M, effective dissociation constant controlling coverage
beta_block = 25;       % dimensionless, how strongly target binding increases Rct
results = table('Size',[numel(Ctarget_labels),9], ...
                  'VariableTypes',{'string','double',...}, ...
                  'VariableNames',{'Label','Ct',...});
  results{k,:} = {string(Ctarget_labels{k}), Ct, ...}; % works if types match
results{k,1} = string(Ctarget_labels{k});
  results{k,2} = Ct;
  % etc.
results = cell(numel(Ctarget_labels), 9); % example: 9 columns
results = cell(numel(Ctarget_labels), 9); % example: 9 columns
results = cell(numel(Ctarget_labels), 9); % example: 9 columns
results = cell(numel(Ctarget_labels), 9); % example: 9 columns
results = cell(numel(Ctarget_labels), 9); % example: 9 columns
gamma_Cdl = 0.30;      % fraction, how much Cdl decreases at saturation (0 to <1)

% Concentrations of target miRNA for visualization/calibration
Ctarget_list = [0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11];   % M
Ctarget_labels = {'Blank', '1 fM', '10 fM', '100 fM', '1 pM', '10 pM'};
colors = lines(numel(Ctarget_list));

% Toggle Warburg tail on/off for comparison
useWarburg = true;

%% ===================== CORE ELECTROCHEMICAL VALUES =====================
% Exchange current and intrinsic Rct from kinetics + redox probe concentration
% i0 = n*F*A*k0*C*
i0 = n*F*A*k0*Credox_m3;      % A

% Butler-Volmer current for the chosen DC overpotential (for information only)
i_BV = i0*( exp((1-alpha)*n*F*eta_dc/(R*T)) - exp(-alpha*n*F*eta_dc/(R*T)) );

% Intrinsic charge-transfer resistance before target binding
Rct0 = (R*T)/(n*F*i0);        % ohm

%% ======================== BUILD IMPEDANCE CURVES ========================
Z_noW = cell(numel(Ctarget_list),1);   % no Warburg -> clean semicircle
Z_withW = cell(numel(Ctarget_list),1); % with Warburg -> tail at low f

results = table('Size',[numel(Ctarget_list), 9], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Label','Ctarget_M','Theta','Rct0_ohm','RctEff_ohm','CdlEff_uF','Tau_s','fPeak_Hz','i0_uA'});

for k = 1:numel(Ctarget_list)
    Ct = Ctarget_list(k);
    theta = Ct/(Ct + Kd_target + eps);                    % simple Langmuir-like surface coverage
    Rct_eff = Rct0*(1 + beta_block*theta);               % target binding increases charge-transfer resistance
    Cdl_eff = Cdl0*(1 - gamma_Cdl*theta);                % target binding decreases effective capacitance slightly
    Cdl_eff = max(Cdl_eff, 1e-9);                        % protect against nonphysical values
    tau = Rct_eff*Cdl_eff;
    fPeak = 1/(2*pi*tau);                                % ideal RC estimate

    Zw = sigmaW*(1 - 1i)./sqrt(w);                       % semi-infinite Warburg

    % Randles circuit without Warburg: Rs + [ Cdl || Rct ]
    Zfar_noW = Rct_eff*ones(size(w));
    Zpar_noW = 1./(1./Zfar_noW + 1i*w*Cdl_eff);
    Z_noW{k} = Rs + Zpar_noW;

    % Randles circuit with Warburg: Rs + [ Cdl || (Rct + Zw) ]
    Zfar_withW = Rct_eff + Zw;
    Zpar_withW = 1./(1./Zfar_withW + 1i*w*Cdl_eff);
    Z_withW{k} = Rs + Zpar_withW;

    results{k, :} = {string(Ctarget_labels{k}), Ct, theta, Rct0, Rct_eff, Cdl_eff*1e6, tau, fPeak, i0*1e6};
end

%% ========================= CONSOLE OUTPUT =========================
fprintf('\n================ EDUCATIONAL EIS / NYQUIST MODEL ================\n');
fprintf('Small AC perturbation amplitude (Vac) = %.3f mV\n', Vac*1e3);
fprintf('Frequency range = %.2e Hz to %.2e Hz\n', fmin, fmax);
fprintf('Redox-probe concentration = %.2f mM\n', Credox_mM);
fprintf('k0 = %.2e m/s\n', k0);
fprintf('Area = %.3e m^2\n', A);
fprintf('Exchange current i0 = %.3f uA\n', i0*1e6);
fprintf('Intrinsic Rct0 = %.2f ohm\n', Rct0);
fprintf('Butler-Volmer current at eta_dc = %.3f V -> i = %.3f uA\n', eta_dc, i_BV*1e6);
fprintf('\nINTERPRETATION:\n');
fprintf('  * Higher k0  -> higher i0  -> smaller Rct -> smaller semicircle\n');
fprintf('  * Higher redox-probe concentration -> higher i0 -> smaller Rct\n');
fprintf('  * Higher target concentration (more blocking) -> larger Rct -> larger semicircle\n');
fprintf('  * Warburg tail appears at LOW frequencies\n');
fprintf('  * Left side of Nyquist = HIGH frequency\n');
fprintf('  * Right side / tail   = LOW frequency\n');
fprintf('===============================================================\n\n');

disp(results);

%% ====================== CHOOSE A CURVE TO ANNOTATE ======================
% Use a mid/high target concentration so changes are visible
idxPlot = min(5, numel(Ctarget_list));     % e.g., 1 pM
Zp_noW = Z_noW{idxPlot};
Zp_W = Z_withW{idxPlot};

% Find representative points for annotations
[~, idxHighF] = min(abs(f - 1e5));
[~, idxMidF]  = min(abs(f - results.fPeak_Hz(idxPlot)));  % approx peak frequency
[~, idxLowF]  = min(abs(f - 0.03));
[~, idxImMax] = max(-imag(Zp_W));

%% =============================== FIGURES ===============================
figure('Color','w','Position',[80 80 1400 850]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% ---- (1) Nyquist without Warburg: clean Rs + Rct semicircle ----
nexttile(1); hold on; box on; grid on;
for k = 1:numel(Ctarget_list)
    Zk = Z_noW{k};
    plot(real(Zk), -imag(Zk), 'LineWidth', 2, 'Color', colors(k,:));
end
xlabel('Z'' (\Omega)');
ylabel('-Z'''' (\Omega)');
title('Nyquist sin Warburg (para leer R_s y R_{ct})');
legend(Ctarget_labels, 'Location', 'best');
axis tight;

% annotate Rs and Rct for chosen curve
Rs_left = real(Zp_noW(1));
Rright = real(Zp_noW(end));
Rct_read = Rright - Rs_left;
plot(Rs_left,0,'ko','MarkerFaceColor','k');
plot(Rright,0,'ko','MarkerFaceColor','k');
text(Rs_left, 0, '  R_s (alta f)', 'VerticalAlignment','bottom', 'FontWeight','bold');
text(Rright, 0, '  R_s + R_{ct}', 'VerticalAlignment','top', 'FontWeight','bold');
annotation('doublearrow',[0.21 0.31],[0.53 0.53]); %#ok<*UNRCH>
text(mean([Rs_left Rright]), max(-imag(Zp_noW))*0.08, sprintf('R_{ct} = %.1f \Omega', Rct_read), ...
    'HorizontalAlignment','center', 'FontWeight','bold', 'BackgroundColor','w');

% ---- (2) Nyquist with Warburg: show low-frequency tail ----
nexttile(2); hold on; box on; grid on;
for k = 1:numel(Ctarget_list)
    Zk = Z_withW{k};
    plot(real(Zk), -imag(Zk), 'LineWidth', 2, 'Color', colors(k,:));
end
xlabel('Z'' (\Omega)');
ylabel('-Z'''' (\Omega)');
title('Nyquist con Warburg (cola difusiva a baja f)');
legend(Ctarget_labels, 'Location', 'best');
axis tight;

% frequency annotations on the chosen curve
plot(real(Zp_W(idxHighF)), -imag(Zp_W(idxHighF)), 'ko', 'MarkerFaceColor','k');
plot(real(Zp_W(idxImMax)), -imag(Zp_W(idxImMax)), 'ko', 'MarkerFaceColor','k');
plot(real(Zp_W(idxLowF)), -imag(Zp_W(idxLowF)), 'ko', 'MarkerFaceColor','k');
text(real(Zp_W(idxHighF)), -imag(Zp_W(idxHighF)), '  Alta f \rightarrow R_s', 'FontWeight','bold');
text(real(Zp_W(idxImMax)), -imag(Zp_W(idxImMax)), '  -Z'''' máximo', 'FontWeight','bold');
text(real(Zp_W(idxLowF)), -imag(Zp_W(idxLowF)), '  Baja f \rightarrow difusión (Warburg)', 'FontWeight','bold');

% ---- (3) Bode magnitude ----
nexttile(3); hold on; box on; grid on;
for k = 1:numel(Ctarget_list)
    Zk = Z_withW{k};
    semilogx(f, abs(Zk), 'LineWidth', 2, 'Color', colors(k,:));
end
set(gca,'XDir','reverse');  % high frequency at left, low at right to match Nyquist reading intuition
xlabel('Frecuencia (Hz)');
ylabel('|Z| (\Omega)');
title('Bode magnitud');
legend(Ctarget_labels, 'Location', 'best');

% ---- (4) Bode phase ----
nexttile(4); hold on; box on; grid on;
for k = 1:numel(Ctarget_list)
    Zk = Z_withW{k};
    semilogx(f, angle(Zk)*180/pi, 'LineWidth', 2, 'Color', colors(k,:));
end
set(gca,'XDir','reverse');
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
title('Bode fase');
legend(Ctarget_labels, 'Location', 'best');

sgtitle('Modelo visual de EIS / Randles para biosensor CoNi-MOF', 'FontWeight','bold');

%% ======================= CALIBRATION CURVE FIGURE =======================
% Educational calibration using Rct increase vs target concentration
blank_Rct = results.RctEff_ohm(1);
DeltaRct = results.RctEff_ohm - blank_Rct;

% Exclude blank for log plot
Ct_nonzero = results.Ctarget_M(2:end);
Delta_nonzero = DeltaRct(2:end);
logC = log10(Ct_nonzero);

% Linear fit in log concentration
p = polyfit(logC, Delta_nonzero, 1);
fitLine = polyval(p, logC);
R2 = 1 - sum((Delta_nonzero - fitLine).^2)/sum((Delta_nonzero - mean(Delta_nonzero)).^2);

figure('Color','w','Position',[180 120 780 560]);
plot(logC, Delta_nonzero, 'o', 'MarkerSize', 8, 'LineWidth', 2); hold on; grid on; box on;
plot(logC, fitLine, '-', 'LineWidth', 2);
for k = 2:numel(Ctarget_list)
    text(logC(k-1), Delta_nonzero(k-1), ['  ' Ctarget_labels{k}], 'FontWeight','bold');
end
xlabel('log_{10}(C_{target} / M)');
ylabel('\DeltaR_{ct} (\Omega)');
title('Curva de calibración educativa: \DeltaR_{ct} vs log C_{target}');
legend('Datos simulados','Ajuste lineal','Location','best');

eqnStr = sprintf('\\DeltaR_{ct} = %.3f log_{10}(C) %+0.3f    (R^2 = %.4f)', p(1), p(2), R2);
text(min(logC)+0.05, max(Delta_nonzero)*0.85, eqnStr, 'BackgroundColor','w', 'EdgeColor','k');

%% ======================== OPTIONAL PARAMETER NOTES ======================
% Try changing these and rerunning:
%   k0          -> larger k0 makes Rct smaller (faster kinetics, smaller semicircle)
%   Credox_mM   -> larger redox concentration makes Rct smaller
%   Rs          -> shifts the whole Nyquist plot to the right
%   Cdl0        -> changes the apex frequency and shape of the arc
%   sigmaW      -> controls how visible the low-frequency Warburg tail is
%   Kd_target   -> controls how quickly target concentration increases coverage
%   beta_block  -> controls how strongly target binding increases Rct
%   fmin        -> if you lower fmin, the Warburg tail becomes more visible
%   fmax        -> if you raise fmax, you see the left high-frequency intercept better

