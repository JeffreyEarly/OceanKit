function summarizeSourcesSinksReservoirs(self,options)
% Summarize Sources Sinks Reservoirs.
%
% summarizeSourcesSinksReservoirs is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Summaries
% - Declaration: summarizeSourcesSinksReservoirs(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Parameter shouldDisplayExactValues: (optional) input argument `shouldDisplayExactValues` (default: false)
% - Parameter shouldDisplayDiffValues: (optional) input argument `shouldDisplayDiffValues` (default: false)
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
    options.shouldDisplayExactValues = false
    options.shouldDisplayDiffValues = false
end

[sources, sinks, inertial_tx, inertial_cascade, ddt, energy] = self.filterEnergyForSourcesSinksReservoirs(timeIndices=options.timeIndices);


%
% Sources table
%
sourcesTable = struct2table(sources);
sourcesTable.Properties.RowNames = sourcesTable.fancyName;
sourcesTable.name = [];
sourcesTable.fancyName = [];
sourcesTable.te_gmda = sourcesTable.te_gmda/self.flux_scale;
sourcesTable.te_wave = sourcesTable.te_wave/self.flux_scale;
sourcesTable.te_damp = sourcesTable.te_damp/self.flux_scale;

if options.shouldDisplayExactValues
    sourcesTable.te_exact = sourcesTable.te_exact/self.flux_scale;
    sourcesTable.te_exact_damp = sourcesTable.te_exact_damp/self.flux_scale;
    if options.shouldDisplayDiffValues
        sourcesTable.te_diff = sourcesTable.te_exact - (sourcesTable.te_gmda + sourcesTable.te_wave);
        sourcesTable.te_damp_diff = sourcesTable.te_exact_damp - sourcesTable.te_damp;
    end
else
    sourcesTable.te_exact=[];
    sourcesTable.te_exact_damp=[];
end
sourcesTable = addTotalsRow(sourcesTable);
disp(sourcesTable)

%
% Sinks table
%
sinksTable = struct2table(sinks);
sinksTable.Properties.RowNames = sinksTable.fancyName;
sinksTable.name = [];
sinksTable.fancyName = [];
sinksTable.te_gmda = sinksTable.te_gmda/self.flux_scale;
sinksTable.te_wave = sinksTable.te_wave/self.flux_scale;
sinksTable.te_damp = sinksTable.te_damp/self.flux_scale;

if options.shouldDisplayExactValues
    sinksTable.te_exact = sinksTable.te_exact/self.flux_scale;
    sinksTable.te_exact_damp = sinksTable.te_exact_damp/self.flux_scale;
    if options.shouldDisplayDiffValues
        sinksTable.te_diff = sinksTable.te_exact - (sinksTable.te_gmda + sinksTable.te_wave);
        sinksTable.te_damp_diff = sinksTable.te_exact_damp - sinksTable.te_damp;
    end
else
    sinksTable.te_exact=[];
    sinksTable.te_exact_damp=[];
end
sinksTable   = addTotalsRow(sinksTable);
disp(sinksTable)

%
% Grand total (sources+sinks) table
%
% Extract the Total rows
srcTotal = sourcesTable("Total", :);
snkTotal = sinksTable("Total", :);

% Ensure same variable order (important if tables were built separately)
snkTotal = snkTotal(:, sourcesTable.Properties.VariableNames);

% Identify numeric variables
numVars = varfun(@isnumeric, sourcesTable, "OutputFormat","uniform");

% Create a new one-row table of the same size
grandTotal = srcTotal;
grandTotal(:,:) = {missing};   % fill everything with missing to start

% Add numeric columns
for v = sourcesTable.Properties.VariableNames(numVars)
    vn = v{1};
    grandTotal.(vn) = srcTotal.(vn) + snkTotal.(vn);
end

% Give the row a clear name
grandTotal.Properties.RowNames = "GrandTotal";

disp(grandTotal)

%
% Inertial table
%
inertialTable = struct2table(inertial_tx);
if ~options.shouldDisplayExactValues
    inertialTable.Properties.RowNames = inertialTable.name;
    inertialTable("te_exact",:) = [];
end
inertialTable.Properties.RowNames = inertialTable.fancyName;
inertialTable.name = [];
inertialTable.fancyName = [];
inertialTable.te_gmda = inertialTable.te_gmda/self.flux_scale;
inertialTable.te_wave = inertialTable.te_wave/self.flux_scale;
inertialTable.te_damp = inertialTable.te_damp/self.flux_scale;

disp(inertialTable)

%
% Inertial table
%
inertialTable = struct2table(inertial_cascade);
if ~options.shouldDisplayExactValues
    inertialTable.Properties.RowNames = inertialTable.name;
    inertialTable("te_exact",:) = [];
end
inertialTable.Properties.RowNames = inertialTable.fancyName;
inertialTable.name = [];
inertialTable.fancyName = [];
inertialTable.te_gmda = inertialTable.te_gmda/self.flux_scale;
inertialTable.te_wave = inertialTable.te_wave/self.flux_scale;
inertialTable.te_damp = inertialTable.te_damp/self.flux_scale;

disp(inertialTable)

%
% ddt table
%
ddtTable = struct2table(ddt);
ddtTable.te_gmda = ddtTable.te_gmda/self.flux_scale;
ddtTable.te_wave = ddtTable.te_wave/self.flux_scale;
ddtTable.te_damp = ddtTable.te_damp/self.flux_scale;
if options.shouldDisplayExactValues
    ddtTable.te_exact = ddtTable.te_exact/self.flux_scale;
    if options.shouldDisplayDiffValues
        ddtTable.te_diff = ddtTable.te_exact - (ddtTable.te_gmda + ddtTable.te_wave);
    end
else
    ddtTable.te_exact=[];
end
disp(ddtTable)



function T = addTotalsRow(T)
    % Build a 1-row table of column sums (works for numeric vars)
    totalRow = varfun(@sum, T, "OutputFormat","table");

    % Rename summed variables back to the originals (remove 'sum_' prefix)
    totalRow.Properties.VariableNames = T.Properties.VariableNames;

    % Give the totals row a unique RowName to avoid duplicates
    totalRow.Properties.RowNames = matlab.lang.makeUniqueStrings( ...
        "Total", T.Properties.RowNames, namelengthmax);

    % Append
    T = [T; totalRow];
end

end