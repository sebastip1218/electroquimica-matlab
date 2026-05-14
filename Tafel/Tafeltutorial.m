clear; clc; close all;

%% ============================================================
%  TUTORIAL VISUAL: efecto de alpha y j0
% =============================================================

% Constantes
n = 1;
F = 96485;
R = 8.314;
T = 298.15;
f = n*F/(R*T);

% Malla de sobrepotencial
eta = linspace(-0.4, 0.4, 1000);

%% ============================================================
%  PARTE A: efecto de alpha
% =============================================================

j0_fixed = 1e-3;                % A/cm^2
alphas = [0.2 0.5 0.8];

figure('Color','w');

subplot(2,2,1); hold on;
for k = 1:length(alphas)
    alpha = alphas(k);
    j_an  =  j0_fixed .* exp((1-alpha).*f.*eta);
    j_cat = -j0_fixed .* exp(-alpha.*f.*eta);
    j_net = j_an + j_cat;
    plot(eta, j_net, 'LineWidth', 2, ...
        'DisplayName', sprintf('\\alpha = %.1f', alpha));
end
xline(0,'k:'); yline(0,'k:');
grid on;
xlabel('\eta / V');
ylabel('j / A cm^{-2}');
title('Butler-Volmer: efecto de \alpha');
legend('Location','best');
set(gca,'FontSize',12,'LineWidth',1.1);

subplot(2,2,2); hold on;
for k = 1:length(alphas)
    alpha = alphas(k);
    j_an  =  j0_fixed .* exp((1-alpha).*f.*eta);
    plot(eta, log(abs(j_an)), 'LineWidth', 2, ...
        'DisplayName', sprintf('anódica, \\alpha = %.1f', alpha));
end
grid on;
xlabel('\eta / V');
ylabel('ln|j_{an}|');
title('Tafel anódica: efecto de \alpha');
legend('Location','best');
set(gca,'FontSize',12,'LineWidth',1.1);

subplot(2,2,3); hold on;
for k = 1:length(alphas)
    alpha = alphas(k);
    j_cat = -j0_fixed .* exp(-alpha.*f.*eta);
    plot(eta, log(abs(j_cat)), 'LineWidth', 2, ...
        'DisplayName', sprintf('catódica, \\alpha = %.1f', alpha));
end
grid on;
xlabel('\eta / V');
ylabel('ln|j_{cat}|');
title('Tafel catódica: efecto de \alpha');
legend('Location','best');
set(gca,'FontSize',12,'LineWidth',1.1);

subplot(2,2,4); hold on;
for k = 1:length(alphas)
    alpha = alphas(k);
    y_an  = log(j0_fixed) + (1-alpha).*f.*eta;
    y_cat = log(j0_fixed) - alpha.*f.*eta;
    plot(eta, y_an, '--', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Tafel anódica, \\alpha = %.1f', alpha));
    plot(eta, y_cat, '-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Tafel catódica, \\alpha = %.1f', alpha));
end
grid on;
xlabel('\eta / V');
ylabel('ln|j|');
title('Rectas de Tafel ideales');
legend('Location','eastoutside');
set(gca,'FontSize',12,'LineWidth',1.1);

sgtitle('Tutorial visual: efecto de \alpha');

%% ============================================================
%  PARTE B: efecto de j0
% =============================================================

alpha_fixed = 0.5;
j0_values = [1e-6 1e-4 1e-2];

figure('Color','w');

subplot(2,2,1); hold on;
for k = 1:length(j0_values)
    j0 = j0_values(k);
    j_an  =  j0 .* exp((1-alpha_fixed).*f.*eta);
    j_cat = -j0 .* exp(-alpha_fixed.*f.*eta);
    j_net = j_an + j_cat;
    plot(eta, j_net, 'LineWidth', 2, ...
        'DisplayName', sprintf('j0 = %.0e', j0));
end
xline(0,'k:'); yline(0,'k:');
grid on;
xlabel('\eta / V');
ylabel('j / A cm^{-2}');
title('Butler-Volmer: efecto de j_0');
legend('Location','best');
set(gca,'FontSize',12,'LineWidth',1.1);

subplot(2,2,2); hold on;
for k = 1:length(j0_values)
    j0 = j0_values(k);
    j_net = j0 .* (exp((1-alpha_fixed).*f.*eta) - exp(-alpha_fixed.*f.*eta));
    plot(eta, log(abs(j_net)+1e-30), 'LineWidth', 2, ...
        'DisplayName', sprintf('j0 = %.0e', j0));
end
grid on;
xlabel('\eta / V');
ylabel('ln|j_{net}|');
title('Curva logarítmica: efecto de j_0');
legend('Location','best');
set(gca,'FontSize',12,'LineWidth',1.1);

subplot(2,2,3); hold on;
for k = 1:length(j0_values)
    j0 = j0_values(k);
    y_an  = log(j0) + (1-alpha_fixed).*f.*eta;
    y_cat = log(j0) - alpha_fixed.*f.*eta;
    plot(eta, y_an, '--', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('anódica, j0 = %.0e', j0));
    plot(eta, y_cat, '-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('catódica, j0 = %.0e', j0));
end
grid on;
xlabel('\eta / V');
ylabel('ln|j|');
title('Rectas de Tafel: efecto de j_0');
legend('Location','eastoutside');
set(gca,'FontSize',12,'LineWidth',1.1);

subplot(2,2,4); hold on;
for k = 1:length(j0_values)
    j0 = j0_values(k);
    intercepto = log(j0);
    plot(k, intercepto, 'o', 'MarkerSize', 10, 'LineWidth', 2, ...
        'DisplayName', sprintf('j0 = %.0e', j0));
    text(k+0.05, intercepto, sprintf('ln(j0)=%.2f', intercepto), 'FontSize', 11);
end
grid on;
xlabel('Caso');
ylabel('ln(j_0)');
title('Intercepto de Tafel');
set(gca,'XTick',1:length(j0_values),'XTickLabel',{'1e-6','1e-4','1e-2'});
set(gca,'FontSize',12,'LineWidth',1.1);

sgtitle('Tutorial visual: efecto de j_0');
