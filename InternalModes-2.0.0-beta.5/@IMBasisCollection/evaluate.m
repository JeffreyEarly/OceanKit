function values = evaluate(self,z,options)
% Evaluate continuous structures with bounded intermediate allocations.
%
% Returns nZ x nColumns x nSelectedPages, preserving page/column order and
% repeated selections. Empty pages/columns selects all. The final output is
% allocated in full; use evaluateChunks to stream without that allocation.
% sampleChunkSize and pageChunkSize bound evaluation temporaries.
%
% - Topic: Evaluate collections
% - Parameter z: physical coordinate column
% - Parameter options.variable: scalar variable name
% - Parameter options.derivativeOrder: supported vertical derivative order
% - Parameter options.pages: requested-page indices, empty means all
% - Parameter options.columns: explicit scientific column indices
% - Parameter options.sampleChunkSize: maximum rows evaluated together
% - Parameter options.pageChunkSize: maximum pages evaluated together
% - Returns values: continuous basis values by requested page
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.pages (1,:) double {mustBeInteger,mustBePositive} = []
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
    options.sampleChunkSize (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 1024
    options.pageChunkSize (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 16
end
arguments (Output)
    values double
end
pages = options.pages;
if isempty(pages)
    pages = 1:numel(self.basisIndex);
elseif any(pages > numel(self.basisIndex))
    error("IMBasisCollection:InvalidPage","pages must select existing requested pages.");
end
nColumns = numel(options.columns);
if isempty(options.columns)
    nColumns = numel(self.metadata(self.basisIndex(pages(1))).columnLabels);
end
values = zeros(numel(z),nColumns,numel(pages));
args = namedargs2cell(options);
self.evaluateChunks(z,@acceptChunk,args{:});

    function acceptChunk(chunk,rows,positions,~)
        values(rows,:,positions) = chunk;
    end
end
