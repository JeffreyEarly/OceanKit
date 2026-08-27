classdef (Sealed) WVFourierStorageLayout
    % Describe how horizontal Fourier storage maps to the WaveVortex grid.
    %
    % WVFourierStorageLayout is developer infrastructure shared by horizontal
    % transform backends. It describes indexing only: the class does not
    % execute an FFT, normalize coefficients, own a backend buffer, or select
    % a transform backend. Model code continues to expose coefficients on the
    % canonical WV grid with shape [Nbatch,Nkl].
    %
    % Four representations appear in this contract:
    %
    % * Spatial arrays have shape [Nx,Ny,Nbatch].
    % * Fourier storage has its natural backend shape. Full-complex storage is
    %   [Nx,Ny,Nbatch], half-x storage is [floor(Nx/2)+1,Ny,Nbatch], and
    %   half-y storage is [Nx,floor(Ny/2)+1,Nbatch].
    % * A Fourier row view reshapes the first two storage dimensions to
    %   [nFourierStorageRows,Nbatch]. No data is reordered by this reshape.
    % * The WV grid has shape [Nbatch,Nkl] and uses the ordering defined by
    %   WVGeometryDoublyPeriodic.
    %
    % For real spatial fields, Fourier coefficients obey
    %
    % $$\hat u(-k,-l)=\overline{\hat u(k,l)}.$$
    %
    % A full-complex layout stores both sides of this relation. A
    % hermitian-half layout stores one side and recovers WV modes with a
    % negative compressed-direction mode by conjugating the stored positive
    % mode. On zero and even-grid Nyquist boundaries, the other horizontal
    % mode may still require explicit Hermitian completion. Modes that are
    % their own conjugates are made exactly real during insertion.
    %
    % Mapping indices are one-based uint64 column vectors. They never contain
    % Nbatch- or Nz-replicated offsets. A performance-sensitive backend may
    % consume these properties directly, instead of calling the convenience
    % mapping methods, when doing so preserves a measured MATLAB indexing and
    % assignment expression.
    %
    % This visible, sealed class is intended for WaveVortex backend developers.
    % It is not an end-user modeling API and is not a supported subclassing
    % point.
    %
    properties (SetAccess=private)
        % Physical horizontal grid shape [Nx,Ny].
        % - Topic: Describe Fourier storage
        % - Developer: true
        horizontalGridSize

        % Natural two-dimensional Fourier storage shape.
        % Full-complex, half-x, and half-y layouts have shapes [Nx,Ny],
        % [floor(Nx/2)+1,Ny], and [Nx,floor(Ny/2)+1], respectively.
        % - Topic: Describe Fourier storage
        % - Developer: true
        fourierStorageSize

        % Number of rows in the two-dimensional Fourier row view.
        % This is prod(fourierStorageSize).
        % - Topic: Describe Fourier storage
        % - Developer: true
        nFourierStorageRows

        % Number of horizontal coefficients in the canonical WV grid.
        % - Topic: Describe Fourier storage
        % - Developer: true
        Nkl

        % Fourier storage representation, "full-complex" or "hermitian-half".
        % - Topic: Describe Fourier storage
        % - Developer: true
        fourierStorageType

        % Compressed horizontal dimension for Hermitian-half storage.
        % The value is [] for full-complex storage, 1 for half-x storage, or
        % 2 for half-y storage.
        % - Topic: Describe Fourier storage
        % - Developer: true
        compressedDimension

        % Measured MATLAB mapping implementation used by the layout.
        % "two-dimensional-rows" denotes indices into a
        % [nFourierStorageRows,Nbatch] row view.
        % - Topic: Describe Fourier storage
        % - Developer: true
        mappingMethod = "two-dimensional-rows"

        % Exact bytes occupied by all one-based uint64 mapping arrays.
        % - Topic: Inspect Fourier storage
        % - Developer: true
        mappingMemoryBytes

        % Fourier rows copied directly to the corresponding WV indices.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        fourierRowsForDirectWVIndices uint64

        % WV-grid indices supplied directly by stored Fourier rows.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        directWVIndices uint64

        % Fourier rows conjugated while recovering the corresponding WV modes.
        % These occur when a requested WV mode lies outside Hermitian-half
        % storage in the compressed direction.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        fourierRowsForConjugatedWVIndices uint64

        % WV-grid indices recovered by conjugating stored Fourier rows.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        conjugatedWVIndices uint64

        % Destination rows filled from Hermitian partners before an inverse FFT.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        hermitianCompletionRows uint64

        % Stored rows whose conjugates fill hermitianCompletionRows.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        hermitianSourceRows uint64

        % WV-grid indices corresponding to hermitianSourceRows.
        % This lets a measured backend write completion rows directly from the
        % WV grid without first reading a Fourier row.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        hermitianSourceWVIndices uint64

        % Fourier rows representing modes equal to their own Hermitian partner.
        % These rows are forced real during WV-to-Fourier insertion.
        % - Topic: Map Fourier storage and WV grid
        % - Developer: true
        selfConjugateFourierRows uint64
    end

    methods
        function self = WVFourierStorageLayout(wvg,fourierStorageType,options)
            % Create a Fourier-storage mapping for a doubly periodic geometry.
            %
            % Full-complex storage requires compressedDimension=[]. A
            % Hermitian-half layout requires compressedDimension=1 (half-x) or
            % compressedDimension=2 (half-y).
            %
            % - Topic: Describe Fourier storage
            % - Declaration: layout = WVFourierStorageLayout(wvg,fourierStorageType,options)
            % - Parameter wvg: WVGeometryDoublyPeriodic defining horizontal modes and WV ordering
            % - Parameter fourierStorageType: "full-complex" or "hermitian-half"
            % - Parameter options.compressedDimension: [], 1, or 2 as required by the storage type
            % - Returns layout: read-only WVFourierStorageLayout
            % - Developer: true
            arguments
                wvg (1,1) WVGeometryDoublyPeriodic
                fourierStorageType (1,1) string {mustBeMember(fourierStorageType,["full-complex","hermitian-half"])}
                options.compressedDimension double = []
            end
            validateCompressedDimension(fourierStorageType,options.compressedDimension);

            self.fourierStorageType = fourierStorageType;
            self.compressedDimension = options.compressedDimension;
            self.horizontalGridSize = [wvg.Nx wvg.Ny];
            self.Nkl = wvg.Nkl;
            if fourierStorageType == "full-complex"
                self.fourierStorageSize = self.horizontalGridSize;
                implicitHermitianDimension = wvg.conjugateDimension;
            else
                self.fourierStorageSize = self.horizontalGridSize;
                self.fourierStorageSize(options.compressedDimension) = floor(self.fourierStorageSize(options.compressedDimension)/2)+1;
                implicitHermitianDimension = options.compressedDimension;
            end
            self.nFourierStorageRows = prod(self.fourierStorageSize);

            kMode = double(wvg.kMode_wv(:));
            lMode = double(wvg.lMode_wv(:));
            [self.fourierRowsForDirectWVIndices,self.directWVIndices,self.fourierRowsForConjugatedWVIndices,self.conjugatedWVIndices,representedK,representedL] = self.wvModeMappings(kMode,lMode);
            [self.hermitianCompletionRows,self.hermitianSourceRows,self.selfConjugateFourierRows] = self.hermitianBoundaryMappings(representedK,representedL,implicitHermitianDimension);
            self.hermitianSourceWVIndices = self.wvIndicesForDirectFourierRows(self.hermitianSourceRows);
            if self.fourierStorageType == "full-complex"
                % Preserve the measured builtin assignment contract. Valid
                % real transforms already make self-conjugate values real.
                self.selfConjugateFourierRows = uint64.empty(0,1);
            end
            self.mappingMemoryBytes = self.calculateMappingMemoryBytes();
        end

        function wvArray = transformFromFourierStorageToWVGrid(self,fourierStorageRows)
            % Map a Fourier row view to canonical WV-grid ordering.
            %
            % FourierStorageRows has shape [nFourierStorageRows,Nbatch]. The
            % returned complex array has shape [Nbatch,Nkl]. Direct modes are
            % gathered without conjugation; modes outside compressed storage
            % are recovered through Hermitian conjugation.
            %
            % - Topic: Map Fourier storage and WV grid
            % - Declaration: wvArray = transformFromFourierStorageToWVGrid(fourierStorageRows)
            % - Parameter fourierStorageRows: Fourier row view [nFourierStorageRows,Nbatch]
            % - Returns wvArray: canonical WV-grid coefficients [Nbatch,Nkl]
            % - Developer: true
            self.validateFourierStorageRows(fourierStorageRows);
            nBatch = size(fourierStorageRows,2);
            if self.fourierStorageType == "full-complex" && isempty(self.conjugatedWVIndices) && isequal(self.directWVIndices,uint64((1:self.Nkl)'))
                wvArray = fourierStorageRows(self.fourierRowsForDirectWVIndices,:).';
                return
            end
            wvArray = complex(zeros(nBatch,self.Nkl));
            wvArray(:,self.directWVIndices) = fourierStorageRows(self.fourierRowsForDirectWVIndices,:).';
            if ~isempty(self.conjugatedWVIndices)
                wvArray(:,self.conjugatedWVIndices) = conj(fourierStorageRows(self.fourierRowsForConjugatedWVIndices,:).');
            end
        end

        function fourierStorageRows = allocateFourierStorage(self,nBatch)
            % Allocate a zeroed complex Fourier row view.
            %
            % The returned array has shape [nFourierStorageRows,Nbatch]. Use
            % reshapeFourierRowsToStorage when an FFT backend needs the natural
            % [NxStorage,NyStorage,Nbatch] shape.
            %
            % - Topic: Manage Fourier storage
            % - Declaration: fourierStorageRows = allocateFourierStorage(nBatch)
            % - Parameter nBatch: positive number of independent transform batches
            % - Returns fourierStorageRows: zeroed complex Fourier row view
            % - Developer: true
            arguments
                self (1,1) WVFourierStorageLayout
                nBatch (1,1) double {mustBeInteger,mustBePositive}
            end
            fourierStorageRows = complex(zeros(self.nFourierStorageRows,nBatch));
        end

        function fourierStorageRows = transformFromWVGridToFourierStorage(self,fourierStorageRows,wvArray)
            % Insert WV-grid coefficients into caller-owned Fourier storage.
            %
            % The first input is a caller-owned row view with shape
            % [nFourierStorageRows,Nbatch]. WVArray has shape [Nbatch,Nkl].
            % Direct and conjugated mappings are inserted, required Hermitian
            % boundary rows are completed, and self-conjugate rows are made
            % real. MATLAB may detach a shared input through copy-on-write, so
            % callers must capture and reassign the returned storage:
            %
            %   rows = layout.transformFromWVGridToFourierStorage(rows,wvArray);
            %
            % - Topic: Map Fourier storage and WV grid
            % - Declaration: fourierStorageRows = transformFromWVGridToFourierStorage(fourierStorageRows,wvArray)
            % - Parameter fourierStorageRows: caller-owned Fourier row view [nFourierStorageRows,Nbatch]
            % - Parameter wvArray: canonical WV-grid coefficients [Nbatch,Nkl]
            % - Returns fourierStorageRows: updated Fourier row view; always reassign this value
            % - Developer: true
            self.validateFourierStorageRows(fourierStorageRows);
            if size(wvArray,1) ~= size(fourierStorageRows,2) || size(wvArray,2) ~= self.Nkl
                error("WaveVortexModel:InvalidWVGridShape","WV-grid coefficients must have shape [Nbatch,Nkl]=[%d,%d].",size(fourierStorageRows,2),self.Nkl);
            end
            fourierStorageRows(self.fourierRowsForDirectWVIndices,:) = wvArray(:,self.directWVIndices).';
            if ~isempty(self.fourierRowsForConjugatedWVIndices)
                fourierStorageRows(self.fourierRowsForConjugatedWVIndices,:) = conj(wvArray(:,self.conjugatedWVIndices).');
            end
            if ~isempty(self.hermitianCompletionRows)
                fourierStorageRows(self.hermitianCompletionRows,:) = conj(fourierStorageRows(self.hermitianSourceRows,:));
            end
            if ~isempty(self.selfConjugateFourierRows)
                fourierStorageRows(self.selfConjugateFourierRows,:) = real(fourierStorageRows(self.selfConjugateFourierRows,:));
            end
        end

        function fourierStorageRows = reshapeFourierStorageToRows(self,fourierStorage)
            % Reshape natural Fourier storage to its two-dimensional row view.
            %
            % FourierStorage must begin with fourierStorageSize. All remaining
            % dimensions are combined as Nbatch. Reshape does not reorder data.
            %
            % - Topic: Manage Fourier storage
            % - Declaration: fourierStorageRows = reshapeFourierStorageToRows(fourierStorage)
            % - Parameter fourierStorage: natural Fourier storage [NxStorage,NyStorage,...]
            % - Returns fourierStorageRows: row view [nFourierStorageRows,Nbatch]
            % - Developer: true
            storageSize = size(fourierStorage);
            if numel(storageSize) < 2 || ~isequal(storageSize(1:2),self.fourierStorageSize)
                error("WaveVortexModel:InvalidFourierStorageShape","Fourier storage must begin with shape [%d,%d].",self.fourierStorageSize(1),self.fourierStorageSize(2));
            end
            if mod(numel(fourierStorage),self.nFourierStorageRows) ~= 0
                error("WaveVortexModel:InvalidFourierStorageShape","Fourier storage size is incompatible with the layout.");
            end
            fourierStorageRows = reshape(fourierStorage,self.nFourierStorageRows,[]);
        end

        function fourierStorage = reshapeFourierRowsToStorage(self,fourierStorageRows)
            % Reshape a row view to natural Fourier-storage dimensions.
            %
            % - Topic: Manage Fourier storage
            % - Declaration: fourierStorage = reshapeFourierRowsToStorage(fourierStorageRows)
            % - Parameter fourierStorageRows: Fourier row view [nFourierStorageRows,Nbatch]
            % - Returns fourierStorage: natural storage [NxStorage,NyStorage,Nbatch]
            % - Developer: true
            self.validateFourierStorageRows(fourierStorageRows);
            fourierStorage = reshape(fourierStorageRows,[self.fourierStorageSize size(fourierStorageRows,2)]);
        end

        function ledger = mappingMemoryUsage(self)
            % Return exact memory usage for each mapping array.
            %
            % Each entry records name, MATLAB class, shape, and bytes. The sum
            % of ledger.bytes equals mappingMemoryBytes.
            %
            % - Topic: Inspect Fourier storage
            % - Declaration: ledger = mappingMemoryUsage()
            % - Returns ledger: structure array describing all uint64 mappings
            % - Developer: true
            names = ["fourierRowsForDirectWVIndices" "directWVIndices" "fourierRowsForConjugatedWVIndices" "conjugatedWVIndices" "hermitianCompletionRows" "hermitianSourceRows" "hermitianSourceWVIndices" "selfConjugateFourierRows"];
            ledger = repmat(struct("name","","class","","shape",[],"bytes",0),1,numel(names));
            for iName = 1:numel(names)
                value = self.(names(iName));
                info = whos("value");
                ledger(iName) = struct("name",names(iName),"class",string(class(value)),"shape",double(size(value)),"bytes",double(info.bytes));
            end
        end

    end

    methods (Access=private)
        function [directRows,directIndices,conjugatedRows,conjugatedIndices,representedK,representedL] = wvModeMappings(self,kMode,lMode)
            wvIndices = uint64((1:numel(kMode))');
            representedK = kMode;
            representedL = lMode;
            conjugated = false(size(kMode));
            if self.fourierStorageType == "hermitian-half"
                compressedModes = kMode;
                if self.compressedDimension == 2
                    compressedModes = lMode;
                end
                conjugated = compressedModes < 0;
                representedK(conjugated) = -representedK(conjugated);
                representedL(conjugated) = -representedL(conjugated);
            end
            rows = self.fourierRowIndicesForModeNumbers(representedK,representedL);
            directRows = rows(~conjugated);
            directIndices = wvIndices(~conjugated);
            conjugatedRows = rows(conjugated);
            conjugatedIndices = wvIndices(conjugated);
        end

        function [destinationRows,sourceRows,selfRows] = hermitianBoundaryMappings(self,kMode,lMode,implicitHermitianDimension)
            boundaryModes = lMode;
            boundaryLength = self.horizontalGridSize(2);
            if implicitHermitianDimension == 1
                boundaryModes = kMode;
                boundaryLength = self.horizontalGridSize(1);
            end
            isBoundary = boundaryModes == 0;
            if mod(boundaryLength,2) == 0
                isBoundary = isBoundary | abs(boundaryModes) == boundaryLength/2;
            end
            sourceRows = self.fourierRowIndicesForModeNumbers(kMode(isBoundary),lMode(isBoundary));
            conjugateRows = self.fourierRowIndicesForModeNumbers(-kMode(isBoundary),-lMode(isBoundary));
            selfConjugate = sourceRows == conjugateRows;
            selfRows = unique(sourceRows(selfConjugate),"stable");

            assignedRows = [self.fourierRowsForDirectWVIndices;self.fourierRowsForConjugatedWVIndices];
            needsCompletion = ~selfConjugate & ~ismember(conjugateRows,assignedRows);
            destinationRows = conjugateRows(needsCompletion);
            sourceRows = sourceRows(needsCompletion);
            [destinationRows,uniqueIndices] = unique(destinationRows,"stable");
            sourceRows = sourceRows(uniqueIndices);
        end

        function rows = fourierRowIndicesForModeNumbers(self,kMode,lMode)
            iK = mod(kMode,self.horizontalGridSize(1))+1;
            iL = mod(lMode,self.horizontalGridSize(2))+1;
            if self.fourierStorageType == "hermitian-half"
                if self.compressedDimension == 1
                    iK = min(iK,self.fourierStorageSize(1));
                else
                    iL = min(iL,self.fourierStorageSize(2));
                end
            end
            rows = uint64(iK+self.fourierStorageSize(1)*(iL-1));
        end

        function bytes = calculateMappingMemoryBytes(self)
            ledger = self.mappingMemoryUsage();
            bytes = sum([ledger.bytes]);
        end

        function indices = wvIndicesForDirectFourierRows(self,rows)
            indices = zeros(size(rows),"uint64");
            for iRow = 1:numel(rows)
                iDirect = find(self.fourierRowsForDirectWVIndices == rows(iRow),1);
                if ~isempty(iDirect)
                    indices(iRow) = self.directWVIndices(iDirect);
                end
            end
        end

        function validateFourierStorageRows(self,fourierStorageRows)
            if ~ismatrix(fourierStorageRows) || size(fourierStorageRows,1) ~= self.nFourierStorageRows
                error("WaveVortexModel:InvalidFourierStorageShape","Fourier-storage rows must have shape [%d,Nbatch].",self.nFourierStorageRows);
            end
        end
    end
end

function validateCompressedDimension(fourierStorageType,compressedDimension)
if fourierStorageType == "full-complex"
    if ~isempty(compressedDimension)
        error("WaveVortexModel:InvalidFullFourierStorageLayout","Full-complex storage requires compressedDimension=[].");
    end
elseif ~isscalar(compressedDimension) || ~ismember(compressedDimension,[1 2])
    error("WaveVortexModel:InvalidHalfFourierStorageLayout","Hermitian-half storage requires compressedDimension 1 or 2.");
end
end
