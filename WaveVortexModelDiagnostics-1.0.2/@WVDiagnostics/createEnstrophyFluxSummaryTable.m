function tableString = createEnstrophyFluxSummaryTable(self,options)
% Create a LaTeX table summarizing enstrophy source/sink fluxes.
%
% Create a LaTeX table summarizing enstrophy source/sink fluxes.
% Generates a LaTeX-formatted table listing enstrophy sources, sinks and totals
% using spatial-temporally averaged enstrophy flux diagnostics. The table
% reports both APV (exact) and quadratic (QGPV) flux measures and formats
% columns so ampersands are aligned for human-readable LaTeX output.
%
% - Topic: Summaries
% - Declaration: tableString = createEnstrophyFluxSummaryTable(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter customNames: (optional) dictionary mapping diagnostic names to custom labels (default: configureDictionary("string","string"))
% - Parameter timeIndices: (optional) time indices to average over (default: Inf -> all times)
% - Returns tableString: string containing the LaTeX table
arguments
    self WVDiagnostics
    options.customNames = configureDictionary("string","string")
    % options.fluxTolerance = 1e-2;
    % options.shouldShowUnits = true;
    options.timeIndices = Inf;
end

enstrophy_fluxes = self.exactEnstrophyFluxesSpatialTemporalAverage(timeIndices=options.timeIndices);
enstrophy_fluxes_qgpv = self.quadraticEnstrophyFluxesSpatialTemporalAverage(timeIndices=options.timeIndices);

[Z_quadratic, t] = self.quadraticEnstrophyOverTime(timeIndices=options.timeIndices);
[Z_apv, ~] = self.exactEnstrophyOverTime(timeIndices=options.timeIndices);


sourceIndices = find([enstrophy_fluxes.Z0] > 0);
sinkIndices = find([enstrophy_fluxes.Z0] < 0);
sources = enstrophy_fluxes(sourceIndices);
[~, idx] = sort([sources.Z0],'descend');
sources = sources(idx);
sinks = enstrophy_fluxes(sinkIndices);
[~, idx] = sort([sinks.Z0],'ascend');
sinks = sinks(idx);


tableString = "\n\\begin{tabular}{rlr}\n";
tableString = tableString + "& forcing & flux APV (QGPV) \\\\ \\hline\\hline \n";
% tableString = tableString + "\\hline\n";
enstrophy_fluxes = sources;
qgpvTotal = 0;
for iFlux=1:length(enstrophy_fluxes)
    if isKey(options.customNames,enstrophy_fluxes(iFlux).name)
        fancyName=options.customNames(enstrophy_fluxes(iFlux).name);
    else
        fancyName=enstrophy_fluxes(iFlux).fancyName;
    end
    qgpvIndex = find(enstrophy_fluxes(iFlux).name == [enstrophy_fluxes_qgpv.name]);
    qgpvTotal = qgpvTotal + enstrophy_fluxes_qgpv(qgpvIndex).Z0;
    tableString = tableString + " & " + fancyName + " & ";
    tableString = tableString + sprintf("%.2f (%.2f)",enstrophy_fluxes(iFlux).Z0/self.z_flux_scale,enstrophy_fluxes_qgpv(qgpvIndex).Z0/self.z_flux_scale)  + " \\\\ \n";
end
tableString = tableString + "\\textbf{sources total} & & \\textbf{";
tableString = tableString + sprintf('%.2f (%.2f)}',sum([enstrophy_fluxes.Z0])/self.z_flux_scale,qgpvTotal/self.z_flux_scale)  + " \\\\ \\hline \n";

enstrophy_fluxes = sinks;
qgpvTotal = 0;
for iFlux=1:length(enstrophy_fluxes)
    if isKey(options.customNames,enstrophy_fluxes(iFlux).name)
        fancyName=options.customNames(enstrophy_fluxes(iFlux).name);
    else
        fancyName=enstrophy_fluxes(iFlux).fancyName;
    end
        qgpvIndex = find(enstrophy_fluxes(iFlux).name == [enstrophy_fluxes_qgpv.name]);
    qgpvTotal = qgpvTotal + enstrophy_fluxes_qgpv(qgpvIndex).Z0;
    tableString = tableString + " & " + fancyName + " & ";
    tableString = tableString + sprintf("%.2f (%.2f)",enstrophy_fluxes(iFlux).Z0/self.z_flux_scale,enstrophy_fluxes_qgpv(qgpvIndex).Z0/self.z_flux_scale)  + " \\\\ \n";
end
tableString = tableString + "\\textbf{sinks total} & & \\textbf{";
tableString = tableString + sprintf("%.2f (%.2f)}",sum([enstrophy_fluxes.Z0])/self.z_flux_scale,qgpvTotal/self.z_flux_scale)  + " \\\\ \\hline ";

tableString = tableString + " \\hline\n";
tableString = tableString + "\\textbf{flux total} & & ";
tableString = tableString + sprintf("%.2f (%.2f)",(sum([sources.Z0])+sum([sinks.Z0]))/self.z_flux_scale,sum([enstrophy_fluxes_qgpv.Z0])/self.z_flux_scale)  + " \\\\ \n";

% tableString = tableString + " \\hline\n";
tableString = tableString + "\\textbf{measured} $\\Delta$ & & ";
tableString = tableString + sprintf("%.2f (%.2f)",(Z_apv(end)-Z_apv(1))/(t(end)-t(1))/self.z_flux_scale,(Z_quadratic(end)-Z_quadratic(1))/(t(end)-t(1))/self.z_flux_scale)  + " \\\\ ";

tableString = tableString + " \\hline\n\\end{tabular} \n";

%% ——— align all ampersands for human‐readable LaTeX code ———
% split the tableString into lines
lines = split(string(tableString), "\n");

% identify the data lines (contain '&' but aren't \hline or \end{tabular})
isData = contains(lines, "&") ...
       & ~startsWith(strtrim(lines), "\hline") ...
       & ~contains(lines, "\end{tabular}");

% extract and tokenize those data lines
dataLines  = lines(isData);
nData      = numel(dataLines);
rowTokens  = cell(nData,1);
for i = 1:nData
    rowTokens{i} = strtrim(split(dataLines(i), "&"));
end

% ensure every row has the same number of columns
nCols = max(cellfun(@numel, rowTokens));
for i = 1:nData
    if numel(rowTokens{i}) < nCols
        rowTokens{i}(nCols) = "";
    end
end

colWidths = zeros(1,nCols);
for c = 1:nCols
    % this returns a string array of size [nData×1]
    colText = cellfun(@(r) r(c), rowTokens, 'UniformOutput', true);
    % now strlength operates on a string array
    colWidths(c) = max(strlength(colText));
end

% rebuild each data line with fixed‐width fields
dataIdx = 1;
for i = 1:numel(lines)
    if isData(i)
        cells  = rowTokens{dataIdx};
        padded = strings(1,nCols);
        for c = 1:nCols
            padded(c) = pad(cells(c), colWidths(c), 'right');
        end
        lines(i) = strjoin(padded, " & ");
        dataIdx = dataIdx + 1;
    end
end

% reassemble back into one big string
tableString = strjoin(lines, newline);

% %% ——— align all ampersands for human‐readable LaTeX code ———
% % split into individual lines
% lines = strsplit(sprintf(tableString), "\n");
% % lines = split(tableString, "\n");
% 
% % find only the "data" lines (those with '&', but not pure \hline or \end{tabular})
% isData = cellfun(@(s) ...
%     contains(s, "&") && ...
%     ~startsWith(strtrim(s), "\hline") && ...
%     ~contains(s, "\end{tabular}"), ...
%     lines);
% 
% % extract and trim each cell
% dataCells = {};
% for ii = find(isData)
%     % split on &, trim whitespace
%     tok = strtrim(strsplit(lines{ii}, "&"));
%     dataCells(end+1,1) = {tok};  %#ok<AGROW>
% end
% 
% % figure out how many columns and their max widths
% nCols = max(cellfun(@numel, dataCells));
% % pad any short rows (just in case)
% for i = 1:numel(dataCells)
%     if numel(dataCells{i}) < nCols
%         dataCells{i}(nCols) = {""};
%     end
% end
% colWidths = zeros(1,nCols);
% for c = 1:nCols
%     % compute the longest cell in column c
%     colWidths(c) = max(cellfun(@strlength, cellfun(@(r) r{c}, dataCells, 'UniformOutput',false)));
% end
% 
% % rebuild each data line with fixed-width fields
% idx = 1;
% for ii = 1:numel(lines)
%     if isData(ii)
%         rowCells = dataCells{idx};
%         parts = strings(1,nCols);
%         for c = 1:nCols
%             % left-justify each cell to its column width
%             parts(c) = sprintf("%-*s", colWidths(c), rowCells{c});
%         end
%         % rejoin with ' & '
%         lines{ii} = convertStringsToChars(strjoin(parts, " & "));
%         idx = idx + 1;
%     end
% end
% 
% % reassemble the tableString
% tableString = join(lines, "\n");
end