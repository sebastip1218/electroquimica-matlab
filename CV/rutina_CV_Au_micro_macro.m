clear; clc; close all;

%% ============================================================
%  RUTINA CV: Au microelectrodo vs Au macroelectrodo
%  ============================================================
%  Objetivo:
%  1) Leer automaticamente los .txt exportados por CHI760D.
%  2) Separar los archivos de Au microelectrodo y Au macroelectrodo.
%  3) Tomar el ultimo ciclo estable para calcular corrientes comparables.
%  4) Mostrar visualmente que tramo uso la rutina para el analisis.
%  5) Incluir un diagnostico inspirado en Li 2023: I = Id + Ik + Ic.
%
%  Nota sobre Li 2023:
%  El paper usa EVLS para separar componentes de corriente con distintas
%  dependencias con la velocidad de barrido:
%      Id ~ v^(1/2), Ik ~ v^0, Ic ~ v^1
%  y combina tripletas de voltamperogramas medidas a v/2, v y 2v.
%  Aqui se incluye como diagnostico educativo, no como reemplazo del
%  analisis principal de meseta/picos.

%% 1) AJUSTES EDITABLES
% Este archivo debe estar en:
% C:\Users\sebas\OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey\Profesional\ELECTROQUIMICA
%
% Por eso la carpeta de datos se define relativa al script:
% ELECTROQUIMICA\080526
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd;
end

dataFolder = fullfile(scriptFolder, '080526');
outputFolder = fullfile(scriptFolder, 'resultados_CV');

% En microelectrodo se toma la meseta catodica a potencial bajo.
plateauWindowV = [0.015 0.050];

saveFigures = true;
saveTables = true;
plotOnlyLastCycle = true;         % true = grafica solo el ultimo ciclo estable

% La figura resumen muestra de donde tomo los datos la rutina.
makeSelectionSummaryFigure = true;

% Guarda una figura diagnostica por cada archivo para revisar cada seleccion.
saveSelectionFiguresPerFile = true;
showPerFileSelectionFigures = false;

% Diagnostico inspirado en Li 2023. Requiere tripletas v/2, v, 2v.
makeEvlsDiagnostic = true;
evlsElectrodeType = "Macro";      % "Macro" se parece mas al caso CV con picos del paper
evlsReferenceRates_Vs = [0.05 0.1]; % EVLS tipo articulo: 0.025/0.05/0.1 y 0.05/0.1/0.2 V/s

% Graficas extra con pocas velocidades para que no se encimen las curvas.
makeReducedVoltammograms = true;
% Se usan las mismas 4 velocidades para micro y macro, asi la comparacion
% visual es directa. Elegi el minimo, dos puntos intermedios separados y el
% maximo comun disponible en ambos conjuntos.
reducedScanRates_Vs = [0.00625 0.05 0.2 0.8];

%% 2) PREPARAR CARPETAS
if ~isfolder(dataFolder)
    error("No existe la carpeta de datos: %s", dataFolder);
end

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

files = dir(fullfile(dataFolder, '*.txt'));
if isempty(files)
    error("No se encontraron archivos .txt en: %s", dataFolder);
end

fprintf("Leyendo %d archivos .txt desde:\n%s\n\n", numel(files), dataFolder);

%% 3) LEER Y ANALIZAR CADA ARCHIVO
% Se preasignan columnas para construir una tabla limpia al final.
nFiles = numel(files);

FileName = strings(nFiles, 1);
ElectrodeType = strings(nFiles, 1);
ScanRate_Vs = nan(nFiles, 1);
NameRate_Vs = nan(nFiles, 1);
NameRateMismatch = false(nFiles, 1);
SegmentHeader = nan(nFiles, 1);
Sensitivity_AV = nan(nFiles, 1);
HeaderLine = nan(nFiles, 1);
NumRows = nan(nFiles, 1);
NumSegmentsDetected = nan(nFiles, 1);
LastCycleStart = nan(nFiles, 1);
LastCycleEnd = nan(nFiles, 1);
AnalysisStart = nan(nFiles, 1);
AnalysisEnd = nan(nFiles, 1);
I_ss_A = nan(nFiles, 1);
I_ss_std_A = nan(nFiles, 1);
NPlateau = nan(nFiles, 1);
ip_cat_A = nan(nFiles, 1);
Ep_cat_V = nan(nFiles, 1);
ip_an_A = nan(nFiles, 1);
Ep_an_V = nan(nFiles, 1);
CurrentCharacteristic_A = nan(nFiles, 1);
Warning = strings(nFiles, 1);

cvData = struct("FileName", {}, "ElectrodeType", {}, "ScanRate_Vs", {}, ...
    "E", {}, "I", {}, "Q", {}, "t", {}, "Segments", {}, ...
    "LastCycleStart", {}, "LastCycleEnd", {}, "AnalysisIdx", {}, ...
    "PlateauIdx", {}, "CatPeakIdx", {}, "AnPeakIdx", {});

validCount = 0;

for k = 1:nFiles
    thisFile = fullfile(files(k).folder, files(k).name);
    FileName(k) = string(files(k).name);

    try
        cv = readChiCvFile(thisFile);

        % Clasificacion simple desde el nombre del archivo.
        lowerName = lower(files(k).name);
        if contains(lowerName, 'micro')
            ElectrodeType(k) = "Micro";
        elseif contains(lowerName, 'macro')
            ElectrodeType(k) = "Macro";
        else
            ElectrodeType(k) = "Desconocido";
            Warning(k) = appendWarning(Warning(k), "No se pudo clasificar como Micro o Macro desde el nombre.");
        end

        % La velocidad real se toma del encabezado del instrumento.
        % Esto es importante porque hay archivos cuyo nombre no coincide
        % exactamente con el Scan Rate guardado por CHI.
        nameRate = extractRateFromFileName(files(k).name);
        NameRate_Vs(k) = nameRate;
        ScanRate_Vs(k) = cv.scanRate;
        SegmentHeader(k) = cv.segmentHeader;
        Sensitivity_AV(k) = cv.sensitivity;
        HeaderLine(k) = cv.headerLine;
        NumRows(k) = numel(cv.E);

        segments = splitCvSegments(cv.E);
        NumSegmentsDetected(k) = size(segments, 1);
        [lastStart, lastEnd] = lastCycleBounds(segments, numel(cv.E));
        LastCycleStart(k) = lastStart;
        LastCycleEnd(k) = lastEnd;

        if ~isnan(nameRate) && ~isnan(cv.scanRate)
            rateTolerance = max(1e-10, 0.01 * abs(cv.scanRate));
            if abs(nameRate - cv.scanRate) > rateTolerance
                NameRateMismatch(k) = true;
                Warning(k) = appendWarning(Warning(k), ...
                    sprintf("Nombre sugiere %.5g V/s; encabezado dice %.5g V/s.", nameRate, cv.scanRate));
            end
        end

        analysisIdx = lastStart:lastEnd;
        plateauIdx = [];
        catPeakIdx = NaN;
        anPeakIdx = NaN;

        % Este IF es el corazon de la rutina:
        % Micro = se busca corriente estacionaria aproximada.
        % Macro = se buscan picos catodico y anodico.
        if ElectrodeType(k) == "Micro"
            [iss, issStd, nPlateau, plateauIdx, analysisIdx] = ...
                estimateMicroPlateau(cv.E, cv.I, segments, plateauWindowV);

            I_ss_A(k) = iss;
            I_ss_std_A(k) = issStd;
            NPlateau(k) = nPlateau;
            CurrentCharacteristic_A(k) = iss;

            if nPlateau == 0
                Warning(k) = appendWarning(Warning(k), "No hubo puntos suficientes en la ventana de meseta.");
            end

        elseif ElectrodeType(k) == "Macro"
            [ipCat, epCat, ipAn, epAn, analysisIdx, catPeakIdx, anPeakIdx] = ...
                estimateMacroPeaks(cv.E, cv.I, segments);

            ip_cat_A(k) = ipCat;
            Ep_cat_V(k) = epCat;
            ip_an_A(k) = ipAn;
            Ep_an_V(k) = epAn;
            CurrentCharacteristic_A(k) = ipCat;
        end

        if ~isempty(analysisIdx)
            AnalysisStart(k) = analysisIdx(1);
            AnalysisEnd(k) = analysisIdx(end);
        end

        validCount = validCount + 1;
        cvData(validCount).FileName = files(k).name;
        cvData(validCount).ElectrodeType = char(ElectrodeType(k));
        cvData(validCount).ScanRate_Vs = cv.scanRate;
        cvData(validCount).E = cv.E;
        cvData(validCount).I = cv.I;
        cvData(validCount).Q = cv.Q;
        cvData(validCount).t = cv.t;
        cvData(validCount).Segments = segments;
        cvData(validCount).LastCycleStart = lastStart;
        cvData(validCount).LastCycleEnd = lastEnd;
        cvData(validCount).AnalysisIdx = analysisIdx;
        cvData(validCount).PlateauIdx = plateauIdx;
        cvData(validCount).CatPeakIdx = catPeakIdx;
        cvData(validCount).AnPeakIdx = anPeakIdx;

    catch ME
        ElectrodeType(k) = "Error";
        Warning(k) = appendWarning(Warning(k), "No se pudo leer: " + string(ME.message));
    end
end

summaryTable = table(FileName, ElectrodeType, ScanRate_Vs, NameRate_Vs, ...
    NameRateMismatch, SegmentHeader, Sensitivity_AV, HeaderLine, NumRows, ...
    NumSegmentsDetected, LastCycleStart, LastCycleEnd, AnalysisStart, ...
    AnalysisEnd, I_ss_A, I_ss_std_A, NPlateau, ip_cat_A, Ep_cat_V, ...
    ip_an_A, Ep_an_V, CurrentCharacteristic_A, Warning);

summaryTable = sortrows(summaryTable, ["ElectrodeType", "ScanRate_Vs", "FileName"]);

%% 4) GRAFICAS PRINCIPALES
figMicro = plotVoltammograms(cvData, "Micro", 1e9, "I / nA", ...
    "Voltamperogramas Au microelectrodo: todas las velocidades, ultimo ciclo", plotOnlyLastCycle);
saveFigureIfNeeded(figMicro, outputFolder, 'figura_1_micro_voltamperogramas.png', saveFigures);

figMacro = plotVoltammograms(cvData, "Macro", 1e6, "I / uA", ...
    "Voltamperogramas Au macroelectrodo: todas las velocidades, ultimo ciclo", plotOnlyLastCycle);
saveFigureIfNeeded(figMacro, outputFolder, 'figura_2_macro_voltamperogramas.png', saveFigures);

if makeReducedVoltammograms
    figMicroReduced = plotVoltammogramsSubset(cvData, "Micro", 1e9, "I / nA", ...
        "Au microelectrodo: 4 velocidades seleccionadas", reducedScanRates_Vs, plotOnlyLastCycle);
    saveFigureIfNeeded(figMicroReduced, outputFolder, 'figura_2b_micro_4_velocidades.png', saveFigures);

    figMacroReduced = plotVoltammogramsSubset(cvData, "Macro", 1e6, "I / uA", ...
        "Au macroelectrodo: 4 velocidades seleccionadas", reducedScanRates_Vs, plotOnlyLastCycle);
    saveFigureIfNeeded(figMacroReduced, outputFolder, 'figura_2c_macro_4_velocidades.png', saveFigures);
end

figCurrent = plotCharacteristicCurrents(summaryTable);
saveFigureIfNeeded(figCurrent, outputFolder, 'figura_3_corriente_vs_velocidad.png', saveFigures);

figScaling = plotScalingTests(summaryTable);
saveFigureIfNeeded(figScaling, outputFolder, 'figura_4_pruebas_de_escala.png', saveFigures);

if makeSelectionSummaryFigure
    figSelection = plotSelectionSummary(cvData, plateauWindowV);
    saveFigureIfNeeded(figSelection, outputFolder, 'figura_5_rutina_visual_seleccion.png', saveFigures);
end

if saveSelectionFiguresPerFile
    saveSelectionDiagnostics(cvData, outputFolder, plateauWindowV, saveFigures, showPerFileSelectionFigures);
end

if makeEvlsDiagnostic
    runEvlsArticleDiagnostics(cvData, evlsElectrodeType, evlsReferenceRates_Vs, outputFolder, saveFigures, saveTables);
end

%% 5) GUARDAR RESULTADOS
if saveTables
    writetable(summaryTable, fullfile(outputFolder, 'resumen_CV.csv'));
    save(fullfile(outputFolder, 'resumen_CV.mat'), ...
        'summaryTable', 'cvData', 'dataFolder', 'plateauWindowV');
end

disp(summaryTable(:, ["FileName", "ElectrodeType", "ScanRate_Vs", ...
    "I_ss_A", "ip_cat_A", "ip_an_A", "Warning"]));

fprintf("\nListo. Resultados guardados en:\n%s\n", outputFolder);

%% ============================================================
%  FUNCIONES LOCALES
%  ============================================================

function cv = readChiCvFile(filePath)
    % Lee un archivo .txt de CHI. No depende de que el encabezado tenga
    % siempre el mismo numero de lineas; busca la linea de columnas.
    rawText = fileread(filePath);
    lines = regexp(rawText, '\r\n|\n|\r', 'split')';

    headerLine = find(contains(lines, 'Potential/V') & contains(lines, 'Current/A'), 1, 'first');
    if isempty(headerLine)
        error("No se encontro la linea de columnas Potential/V, Current/A.");
    end

    dataText = strjoin(lines(headerLine + 1:end), newline);
    parsed = textscan(dataText, '%f%f%f%f', ...
        'Delimiter', ',', 'MultipleDelimsAsOne', true, 'CollectOutput', true);

    data = parsed{1};
    if isempty(data) || size(data, 2) < 4
        error("No se encontraron datos numericos de cuatro columnas.");
    end

    cv.E = data(:, 1);
    cv.I = data(:, 2);
    cv.Q = data(:, 3);
    cv.t = data(:, 4);
    cv.headerLine = headerLine;
    cv.scanRate = getHeaderNumber(lines, 'Scan Rate (V/s)');
    cv.segmentHeader = getHeaderNumber(lines, 'Segment');
    cv.sensitivity = getHeaderNumber(lines, 'Sensitivity (A/V)');
end

function value = getHeaderNumber(lines, label)
    % Extrae numeros del encabezado, por ejemplo:
    % Scan Rate (V/s) = 0.1
    value = NaN;
    idx = find(contains(lines, label), 1, 'first');
    if isempty(idx)
        return;
    end

    token = regexp(lines{idx}, '=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', 'tokens', 'once');
    if ~isempty(token)
        value = str2double(token{1});
    end
end

function value = extractRateFromFileName(fileName)
    % Solo sirve como diagnostico. El analisis usa el Scan Rate del encabezado.
    token = regexp(fileName, '\d+(?:\.\d+)?', 'match', 'once');
    if isempty(token)
        value = NaN;
    else
        value = str2double(token);
    end
end

function segments = splitCvSegments(E)
    % Divide el CV cada vez que cambia la direccion del potencial:
    % baja, sube, baja, sube...
    n = numel(E);
    if n < 2
        segments = [1 n];
        return;
    end

    dE = diff(E);
    tol = max(1e-12, 1e-9 * max(abs(E)));

    direction = zeros(size(dE));
    direction(dE > tol) = 1;
    direction(dE < -tol) = -1;

    firstNonZero = find(direction ~= 0, 1, 'first');
    if isempty(firstNonZero)
        segments = [1 n];
        return;
    end

    direction(1:firstNonZero - 1) = direction(firstNonZero);
    for i = firstNonZero + 1:numel(direction)
        if direction(i) == 0
            direction(i) = direction(i - 1);
        end
    end

    turnPoints = find(direction(1:end - 1) .* direction(2:end) < 0) + 1;
    starts = [1; turnPoints(:)];
    ends = [turnPoints(:); n];
    segments = [starts ends];

    minPoints = 5;
    keep = (segments(:, 2) - segments(:, 1) + 1) >= minPoints;
    segments = segments(keep, :);

    if isempty(segments)
        segments = [1 n];
    end
end

function [startIdx, endIdx] = lastCycleBounds(segments, nPoints)
    % El ultimo ciclo estable se toma como los dos ultimos segmentos
    % completos: una rama de ida y una rama de regreso.
    if isempty(segments)
        startIdx = 1;
        endIdx = nPoints;
    elseif size(segments, 1) >= 2
        startIdx = segments(end - 1, 1);
        endIdx = segments(end, 2);
    else
        startIdx = segments(1, 1);
        endIdx = segments(1, 2);
    end
end

function [idx, segmentNumber] = findLastSegmentByDirection(E, segments, desiredDirection, maxSegment)
    % desiredDirection = -1 busca barrido hacia potencial bajo.
    % desiredDirection =  1 busca barrido hacia potencial alto.
    idx = [];
    segmentNumber = NaN;

    if nargin < 4 || isnan(maxSegment)
        maxSegment = size(segments, 1);
    end
    maxSegment = min(maxSegment, size(segments, 1));

    for s = maxSegment:-1:1
        thisIdx = segments(s, 1):segments(s, 2);
        deltaE = E(thisIdx(end)) - E(thisIdx(1));
        if desiredDirection * deltaE > 0
            idx = thisIdx;
            segmentNumber = s;
            return;
        end
    end
end

function [iss, issStd, nPlateau, plateauIdx, analysisIdx] = estimateMicroPlateau(E, I, segments, plateauWindowV)
    % Para microelectrodo: toma la ultima rama catodica y promedia los
    % puntos dentro de la ventana de meseta. Esto aproxima i_ss.
    iss = NaN;
    issStd = NaN;
    nPlateau = 0;
    plateauIdx = [];

    [analysisIdx, ~] = findLastSegmentByDirection(E, segments, -1, size(segments, 1));
    if isempty(analysisIdx)
        [startIdx, endIdx] = lastCycleBounds(segments, numel(E));
        analysisIdx = startIdx:endIdx;
    end

    vMin = min(plateauWindowV);
    vMax = max(plateauWindowV);
    localMask = E(analysisIdx) >= vMin & E(analysisIdx) <= vMax;
    plateauIdx = analysisIdx(localMask);
    values = I(plateauIdx);
    values = values(~isnan(values));

    nPlateau = numel(values);
    if nPlateau > 0
        iss = mean(values);
        issStd = std(values);
    end
end

function [ipCat, epCat, ipAn, epAn, analysisIdx, catPeakIdx, anPeakIdx] = estimateMacroPeaks(E, I, segments)
    % Para macroelectrodo: busca el minimo de la ultima rama catodica y el
    % maximo de la rama anodica inmediatamente anterior.
    ipCat = NaN;
    epCat = NaN;
    ipAn = NaN;
    epAn = NaN;
    catPeakIdx = NaN;
    anPeakIdx = NaN;

    [catIdx, catSegment] = findLastSegmentByDirection(E, segments, -1, size(segments, 1));
    [anIdx, ~] = findLastSegmentByDirection(E, segments, 1, catSegment - 1);

    if isempty(catIdx)
        [startIdx, endIdx] = lastCycleBounds(segments, numel(E));
        catIdx = startIdx:endIdx;
    end
    if isempty(anIdx)
        [startIdx, endIdx] = lastCycleBounds(segments, numel(E));
        anIdx = startIdx:endIdx;
    end

    analysisIdx = min([catIdx(:); anIdx(:)]):max([catIdx(:); anIdx(:)]);

    [ipCat, localCat] = min(I(catIdx));
    catPeakIdx = catIdx(localCat);
    epCat = E(catPeakIdx);

    [ipAn, localAn] = max(I(anIdx));
    anPeakIdx = anIdx(localAn);
    epAn = E(anPeakIdx);
end

function fig = plotVoltammograms(cvData, electrodeType, currentScale, yLabelText, titleText, plotOnlyLastCycle)
    fig = figure('Color', 'w', 'Name', titleText);
    hold on;

    if isempty(cvData)
        title(titleText);
        text(0.5, 0.5, "No hay datos validos.", 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
        return;
    end

    types = string({cvData.ElectrodeType});
    idx = find(types == electrodeType);
    if isempty(idx)
        title(titleText);
        text(0.5, 0.5, "No hay archivos " + electrodeType + ".", 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
        grid on;
        return;
    end

    [~, order] = sort([cvData(idx).ScanRate_Vs]);
    idx = idx(order);
    colors = lines(numel(idx));

    for j = 1:numel(idx)
        item = cvData(idx(j));
        if plotOnlyLastCycle
            useIdx = item.LastCycleStart:item.LastCycleEnd;
        else
            useIdx = 1:numel(item.E);
        end

        label = sprintf('%.5g V/s - %s', item.ScanRate_Vs, erase(string(item.FileName), ".txt"));
        plot(item.E(useIdx), item.I(useIdx) * currentScale, ...
            'LineWidth', 1.3, 'Color', colors(j, :), 'DisplayName', label);
    end

    grid on;
    box on;
    xlabel('E / V');
    ylabel(yLabelText);
    title(titleText);
    legend('Location', 'bestoutside', 'Interpreter', 'none');
    set(gca, 'FontSize', 11, 'LineWidth', 1.0);
end

function fig = plotVoltammogramsSubset(cvData, electrodeType, currentScale, yLabelText, titleText, scanRatesToPlot, plotOnlyLastCycle)
    % Grafica solo unas cuantas velocidades. Esto ayuda a ver la forma de
    % las curvas sin que todos los archivos queden encimados.
    fig = figure('Color', 'w', 'Name', titleText);
    hold on;

    items = filterItemsByType(cvData, electrodeType);
    if isempty(items)
        title(titleText);
        text(0.5, 0.5, "No hay archivos " + electrodeType + ".", ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
        grid on;
        return;
    end

    selected = chooseItemsAtRates(items, scanRatesToPlot);
    colors = lines(numel(selected));

    for j = 1:numel(selected)
        item = selected(j);
        if plotOnlyLastCycle
            useIdx = item.LastCycleStart:item.LastCycleEnd;
        else
            useIdx = 1:numel(item.E);
        end

        label = sprintf('%.5g V/s - %s', item.ScanRate_Vs, erase(string(item.FileName), ".txt"));
        plot(item.E(useIdx), item.I(useIdx) * currentScale, ...
            'LineWidth', 1.7, 'Color', colors(j, :), 'DisplayName', label);
    end

    grid on;
    box on;
    xlabel('E / V');
    ylabel(yLabelText);
    title(titleText);
    legend('Location', 'bestoutside', 'Interpreter', 'none');
    set(gca, 'FontSize', 11, 'LineWidth', 1.0);
end

function fig = plotCharacteristicCurrents(summaryTable)
    fig = figure('Color', 'w', 'Name', 'Corriente caracteristica vs velocidad');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    micro = summaryTable(summaryTable.ElectrodeType == "Micro" & ~isnan(summaryTable.I_ss_A), :);
    macro = summaryTable(summaryTable.ElectrodeType == "Macro" & ~isnan(summaryTable.ip_cat_A), :);

    nexttile;
    if ~isempty(micro)
        micro = sortrows(micro, 'ScanRate_Vs');
        plot(micro.ScanRate_Vs, abs(micro.I_ss_A) * 1e9, 'o-', ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    else
        text(0.5, 0.5, 'Sin datos micro', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    grid on; box on;
    xlabel('v / V s^{-1}');
    ylabel('|i_{ss}| / nA');
    title('Au microelectrodo: meseta catodica');

    nexttile;
    if ~isempty(macro)
        macro = sortrows(macro, 'ScanRate_Vs');
        plot(macro.ScanRate_Vs, abs(macro.ip_cat_A) * 1e6, 's-', ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    else
        text(0.5, 0.5, 'Sin datos macro', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    grid on; box on;
    xlabel('v / V s^{-1}');
    ylabel('|i_{p,cat}| / uA');
    title('Au macroelectrodo: pico catodico');

    sgtitle('Comparacion de corriente caracteristica por tipo de electrodo');
end

function fig = plotScalingTests(summaryTable)
    fig = figure('Color', 'w', 'Name', 'Pruebas de escala');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    micro = summaryTable(summaryTable.ElectrodeType == "Micro" & ~isnan(summaryTable.I_ss_A), :);
    macro = summaryTable(summaryTable.ElectrodeType == "Macro" & ~isnan(summaryTable.ip_cat_A), :);

    nexttile;
    if ~isempty(micro)
        micro = sortrows(micro, 'ScanRate_Vs');
        x = micro.ScanRate_Vs;
        y = abs(micro.I_ss_A) * 1e9;
        plot(x, y, 'o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Datos');
        hold on;
        addLinearFit(x, y, 'Ajuste lineal');
        yline(mean(y), ':', 'Promedio', 'LineWidth', 1.2, 'DisplayName', 'Promedio');
    else
        text(0.5, 0.5, 'Sin datos micro', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    grid on; box on;
    xlabel('v / V s^{-1}');
    ylabel('|i_{ss}| / nA');
    title('Au microelectrodo: |i_{ss}| vs v');
    legend('Location', 'best');

    nexttile;
    if ~isempty(macro)
        macro = sortrows(macro, 'ScanRate_Vs');
        x = sqrt(macro.ScanRate_Vs);
        y = abs(macro.ip_cat_A) * 1e6;
        plot(x, y, 's', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Datos');
        hold on;
        addLinearFit(x, y, 'Ajuste lineal');
    else
        text(0.5, 0.5, 'Sin datos macro', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    grid on; box on;
    xlabel('sqrt(v) / (V s^{-1})^{1/2}');
    ylabel('|i_{p,cat}| / uA');
    title('Au macroelectrodo: |i_p| vs sqrt(v)');
    legend('Location', 'best');

    sgtitle('Pruebas de escala: microelectrodo y macroelectrodo');
end

function fig = plotSelectionSummary(cvData, plateauWindowV)
    % Figura visible para verificar "de donde tomo los datos" la rutina.
    fig = figure('Color', 'w', 'Name', 'Rutina visual de seleccion');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    microIdx = findRepresentativeIndex(cvData, "Micro", 0.1);
    macroIdx = findRepresentativeIndex(cvData, "Macro", 0.1);

    nexttile;
    if ~isnan(microIdx)
        plotSelectionForItem(cvData(microIdx), plateauWindowV, true);
    else
        text(0.5, 0.5, 'Sin ejemplo micro', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end

    nexttile;
    if ~isnan(macroIdx)
        plotSelectionForItem(cvData(macroIdx), plateauWindowV, true);
    else
        text(0.5, 0.5, 'Sin ejemplo macro', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
end

function idx = findRepresentativeIndex(cvData, electrodeType, preferredRate)
    idx = NaN;
    if isempty(cvData)
        return;
    end

    types = string({cvData.ElectrodeType});
    candidates = find(types == electrodeType);
    if isempty(candidates)
        return;
    end

    rates = [cvData(candidates).ScanRate_Vs];
    [~, localIdx] = min(abs(rates - preferredRate));
    idx = candidates(localIdx);
end

function saveSelectionDiagnostics(cvData, outputFolder, plateauWindowV, saveFigures, showFigures)
    if ~saveFigures
        return;
    end

    selectionFolder = fullfile(outputFolder, 'seleccion_por_archivo');
    if ~isfolder(selectionFolder)
        mkdir(selectionFolder);
    end

    for k = 1:numel(cvData)
        if showFigures
            fig = figure('Color', 'w', 'Name', ['Seleccion - ' cvData(k).FileName]);
        else
            fig = figure('Color', 'w', 'Visible', 'off', 'Name', ['Seleccion - ' cvData(k).FileName]);
        end

        plotSelectionForItem(cvData(k), plateauWindowV, false);
        safeName = regexprep(erase(string(cvData(k).FileName), ".txt"), '[^\w.-]', '_');
        saveFigureIfNeeded(fig, selectionFolder, char("seleccion_" + safeName + ".png"), true);

        if ~showFigures
            close(fig);
        end
    end
end

function plotSelectionForItem(item, plateauWindowV, addTitlePrefix)
    % Grafica completa en gris y resalta el tramo usado para analisis.
    if nargin < 3
        addTitlePrefix = false;
    end

    if strcmpi(item.ElectrodeType, 'Micro')
        scale = 1e9;
        yText = 'I / nA';
    else
        scale = 1e6;
        yText = 'I / uA';
    end

    plot(item.E, item.I * scale, '-', 'Color', [0.70 0.70 0.70], ...
        'LineWidth', 0.9, 'DisplayName', 'CV completo');
    hold on;

    if ~isempty(item.AnalysisIdx)
        plot(item.E(item.AnalysisIdx), item.I(item.AnalysisIdx) * scale, ...
            'b-', 'LineWidth', 1.8, 'DisplayName', 'Tramo usado');
    end

    if strcmpi(item.ElectrodeType, 'Micro')
        if ~isempty(item.PlateauIdx)
            scatter(item.E(item.PlateauIdx), item.I(item.PlateauIdx) * scale, ...
                24, 'r', 'filled', 'DisplayName', 'Puntos de meseta');
            yMean = mean(item.I(item.PlateauIdx) * scale);
            yline(yMean, 'r--', 'DisplayName', 'i_{ss} promedio');
        end
        xline(plateauWindowV(1), 'k:', 'DisplayName', 'Ventana meseta');
        xline(plateauWindowV(2), 'k:', 'HandleVisibility', 'off');
    else
        if ~isnan(item.CatPeakIdx)
            scatter(item.E(item.CatPeakIdx), item.I(item.CatPeakIdx) * scale, ...
                55, 'v', 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Pico catodico');
        end
        if ~isnan(item.AnPeakIdx)
            scatter(item.E(item.AnPeakIdx), item.I(item.AnPeakIdx) * scale, ...
                55, '^', 'filled', 'MarkerFaceColor', [0 0.55 0], 'DisplayName', 'Pico anodico');
        end
    end

    grid on; box on;
    xlabel('E / V');
    ylabel(yText);

    if addTitlePrefix
        title(sprintf('%s: %.5g V/s - %s', electrodeLabel(item.ElectrodeType), item.ScanRate_Vs, item.FileName), ...
            'Interpreter', 'none');
    else
        title(sprintf('Seleccion rutina %s: %.5g V/s - %s', electrodeLabel(item.ElectrodeType), item.ScanRate_Vs, item.FileName), ...
            'Interpreter', 'none');
    end

    legend('Location', 'best', 'Interpreter', 'none');
end

function label = electrodeLabel(electrodeType)
    if strcmpi(char(electrodeType), 'Micro')
        label = 'Au microelectrodo';
    elseif strcmpi(char(electrodeType), 'Macro')
        label = 'Au macroelectrodo';
    else
        label = char(electrodeType);
    end
end

function runEvlsArticleDiagnostics(cvData, electrodeType, referenceRates, outputFolder, saveFigures, saveTables)
    % EVLS tipo articulo, basado en EVLS_Bueno.m y Li 2023.
    %
    % La idea teorica es I = Id + Ik + Ic, donde las componentes dependen
    % diferente de la velocidad de barrido:
    %   Id ~ v^(1/2), Ik ~ v^0, Ic ~ v^1
    %
    % Para cada referencia v_ref, se necesita la tripleta:
    %   i_half   = corriente a v_ref/2
    %   i_ref    = corriente a v_ref
    %   i_double = corriente a 2*v_ref
    %
    % Visualmente, como en EVLS_Bueno, se grafican:
    %   i_ref, i_dif = E1, i_kin = E5, i_cap/ads = E6.

    items = filterItemsByType(cvData, electrodeType);
    if isempty(items)
        warning('No hay archivos %s para EVLS.', electrodeType);
        return;
    end

    evlsFolder = fullfile(outputFolder, 'EVLS_tipo_articulo');
    if ~isfolder(evlsFolder)
        mkdir(evlsFolder);
    end

    rates = unique(round([items.ScanRate_Vs], 10));
    availableRefs = findTripletReferenceRates(rates);
    if isempty(availableRefs)
        warning('No hay tripletas v/2, v, 2v disponibles para EVLS.');
        return;
    end

    if isempty(referenceRates) || all(isnan(referenceRates))
        referenceRates = availableRefs;
    end

    usedAny = false;
    for k = 1:numel(referenceRates)
        targetRef = referenceRates(k);
        [distance, closest] = min(abs(availableRefs - targetRef));
        tolerance = max(1e-10, 0.01 * abs(targetRef));
        if distance > tolerance
            warning('No existe tripleta EVLS para v_ref = %.5g V/s. Se omite.', targetRef);
            continue;
        end

        vref = availableRefs(closest);
        halfItem = chooseItemAtRate(items, vref / 2);
        refItem = chooseItemAtRate(items, vref);
        doubleItem = chooseItemAtRate(items, 2 * vref);

        [E_common, i_half, i_ref, i_double] = buildEvlsCommonAxis(halfItem, refItem, doubleItem, 'anodica');
        [i_dif, i_kin, i_capads, E1, E2, E3, E4, E5, E6] = calculateEvlsArticleTerms(i_half, i_ref, i_double);

        fig = plotEvlsArticleFigure(E_common, i_ref, i_dif, i_kin, i_capads, electrodeType, vref);
        fileTag = sprintf('%03dmVs', round(vref * 1000));
        saveFigureIfNeeded(fig, evlsFolder, sprintf('EVLS_tipo_articulo_%s_%s.png', char(electrodeType), fileTag), saveFigures);
        saveFigureIfNeeded(fig, outputFolder, sprintf('figura_%d_evls_tipo_articulo_%s_%s.png', 5 + k, char(electrodeType), fileTag), saveFigures);

        if saveTables
            evlsTable = table(E_common, i_half, i_ref, i_double, i_dif, i_kin, i_capads, E1, E2, E3, E4, E5, E6, ...
                'VariableNames', {'Potential_V','i_half_A','i_ref_A','i_double_A','i_dif_E1_A','i_kin_E5_A','i_capads_E6_A','E1_A','E2_A','E3_A','E4_A','E5_A','E6_A'});
            writetable(evlsTable, fullfile(evlsFolder, sprintf('EVLS_tipo_articulo_%s_%s.csv', char(electrodeType), fileTag)));
        end

        usedAny = true;
    end

    if ~usedAny
        warning('No se genero ninguna grafica EVLS tipo articulo. Revisa las velocidades disponibles.');
    end
end

function vRefs = findTripletReferenceRates(rates)
    vRefs = [];
    for r = rates
        hasHalf = any(abs(rates - r/2) < 1e-10);
        hasDouble = any(abs(rates - 2*r) < 1e-10);
        if hasHalf && hasDouble
            vRefs(end + 1) = r; %#ok<AGROW>
        end
    end
end

function [E_common, i_half, i_ref, i_double] = buildEvlsCommonAxis(halfItem, refItem, doubleItem, branchName)
    % Construye el eje comun de potencial como en EVLS_Bueno:
    % se toma la zona de potencial compartida por las tres velocidades y
    % luego se interpolan las corrientes al mismo eje.
    [E_half, I_half_raw] = getEvlsBranch(halfItem, branchName);
    [E_ref, I_ref_raw] = getEvlsBranch(refItem, branchName);
    [E_double, I_double_raw] = getEvlsBranch(doubleItem, branchName);

    Emin = max([min(E_half), min(E_ref), min(E_double)]);
    Emax = min([max(E_half), max(E_ref), max(E_double)]);
    E_common = linspace(Emin, Emax, 1000)';

    i_half = interp1(E_half, I_half_raw, E_common, 'linear');
    i_ref = interp1(E_ref, I_ref_raw, E_common, 'linear');
    i_double = interp1(E_double, I_double_raw, E_common, 'linear');
end

function [Ebranch, Ibranch] = getEvlsBranch(item, branchName)
    % EVLS_Bueno usa por defecto la rama anodica. Aqui se toma esa misma
    % idea, pero detectando segmentos reales del CV en lugar de asumir un
    % numero fijo de puntos por ciclo.
    if strcmpi(branchName, 'anodica')
        [idx, ~] = findLastSegmentByDirection(item.E, item.Segments, 1, size(item.Segments, 1));
    else
        [idx, ~] = findLastSegmentByDirection(item.E, item.Segments, -1, size(item.Segments, 1));
    end

    if isempty(idx)
        idx = item.LastCycleStart:item.LastCycleEnd;
    end

    Ebranch = item.E(idx);
    Ibranch = item.I(idx);

    [Ebranch, order] = sort(Ebranch);
    Ibranch = Ibranch(order);
    [Ebranch, uniqueIdx] = unique(Ebranch, 'stable');
    Ibranch = Ibranch(uniqueIdx);
end

function [i_dif, i_kin, i_capads, E1, E2, E3, E4, E5, E6] = calculateEvlsArticleTerms(i_half, i_ref, i_double)
    % Misma estructura que EVLS_Bueno.m.
    E1 = -3.4142 .* i_half + 3.4142 .* i_ref;
    E2 =  4.8284 .* i_half - 2.4142 .* i_ref;
    E3 =  3.4142 .* i_half - 2.4142 .* i_ref;

    E4 = -11.6570 .* i_half + 17.4850 .* i_ref - 5.8284 .* i_double;
    E5 =   6.8284 .* i_half -  8.2426 .* i_ref + 2.4142 .* i_double;
    E6 =   4.8284 .* i_half -  8.2426 .* i_ref + 3.4142 .* i_double;

    % Asignacion usada para la grafica tipo articulo solicitada.
    i_dif = E1;
    i_kin = E5;
    i_capads = E6;
end

function fig = plotEvlsArticleFigure(E, i_ref, i_dif, i_kin, i_capads, electrodeType, vref)
    fig = figure('Color', 'w', 'Name', sprintf('EVLS tipo articulo %.5g V/s', vref));
    hold on;

    scale = chooseCurrentScale(electrodeType);
    plot(E, i_ref .* scale, 'LineWidth', 2.2, 'DisplayName', sprintf('CV ref %.0f mV/s', vref * 1000));
    plot(E, i_dif .* scale, 'LineWidth', 2.2, 'DisplayName', 'i dif (E1)');
    plot(E, i_kin .* scale, 'LineWidth', 2.2, 'DisplayName', 'i kin (E5)');
    plot(E, i_capads .* scale, 'LineWidth', 2.2, 'DisplayName', 'i cap/ads (E6)');

    xlabel('E / V', 'FontSize', 14);
    ylabel(currentScaleLabel(electrodeType), 'FontSize', 14);
    title(sprintf('EVLS tipo articulo - %s - referencia %.0f mV/s', electrodeLabel(electrodeType), vref * 1000), 'FontSize', 13);
    legend('Location', 'best');
    set(gca, 'FontSize', 12, 'LineWidth', 1.4);
    box on;
    grid off;
    xlim([min(E), max(E)]);
end
function items = filterItemsByType(cvData, electrodeType)
    items = struct([]);
    if isempty(cvData)
        return;
    end

    types = string({cvData.ElectrodeType});
    idx = find(types == electrodeType);
    if isempty(idx)
        return;
    end

    items = cvData(idx);
end

function item = chooseItemAtRate(items, targetRate)
    rates = [items.ScanRate_Vs];
    [~, idx] = min(abs(rates - targetRate));
    item = items(idx);
end

function selected = chooseItemsAtRates(items, targetRates)
    % Elige un archivo representativo por cada velocidad solicitada.
    % Si hay replicas con la misma velocidad, toma la primera en orden alfabetico.
    if isempty(items)
        selected = items;
        return;
    end

    [~, nameOrder] = sort(string({items.FileName}));
    items = items(nameOrder);
    rates = [items.ScanRate_Vs];
    selected = items([]);

    for k = 1:numel(targetRates)
        [distance, idx] = min(abs(rates - targetRates(k)));
        tolerance = max(1e-10, 0.01 * abs(targetRates(k)));
        if distance <= tolerance
            selected(end + 1) = items(idx); %#ok<AGROW>
        end
    end

    if isempty(selected)
        [~, order] = sort(rates);
        maxItems = min(4, numel(order));
        selected = items(order(1:maxItems));
    end
end

function [Ebranch, Ibranch] = getCathodicBranch(item)
    [idx, ~] = findLastSegmentByDirection(item.E, item.Segments, -1, size(item.Segments, 1));
    if isempty(idx)
        idx = item.LastCycleStart:item.LastCycleEnd;
    end

    Ebranch = item.E(idx);
    Ibranch = item.I(idx);

    % interp1 necesita un eje monotono. Se ordena por potencial.
    [Ebranch, order] = sort(Ebranch);
    Ibranch = Ibranch(order);
    [Ebranch, uniqueIdx] = unique(Ebranch, 'stable');
    Ibranch = Ibranch(uniqueIdx);
end

function Iquery = interpCurrent(Esource, Isource, Equery)
    Iquery = nan(size(Equery));
    validRange = Equery >= min(Esource) & Equery <= max(Esource);
    Iquery(validRange) = interp1(Esource, Isource, Equery(validRange), 'linear');
end

function scale = chooseCurrentScale(electrodeType)
    if electrodeType == "Micro"
        scale = 1e9;
    else
        scale = 1e6;
    end
end

function label = currentScaleLabel(electrodeType)
    if electrodeType == "Micro"
        label = 'I / nA';
    else
        label = 'I / uA';
    end
end

function addLinearFit(x, y, labelText)
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        return;
    end

    p = polyfit(x, y, 1);
    xFit = linspace(min(x), max(x), 100);
    yFit = polyval(p, xFit);
    plot(xFit, yFit, '--', 'LineWidth', 1.3, 'DisplayName', labelText);
end

function saveFigureIfNeeded(fig, outputPath, fileName, saveFigures)
    if ~saveFigures || isempty(fig) || ~isvalid(fig)
        return;
    end

    filePath = fullfile(outputPath, fileName);
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, filePath, 'Resolution', 300);
    else
        saveas(fig, filePath);
    end
end

function out = appendWarning(existing, message)
    existing = string(existing);
    message = string(message);

    if strlength(existing) == 0
        out = message;
    else
        out = existing + " " + message;
    end
end
