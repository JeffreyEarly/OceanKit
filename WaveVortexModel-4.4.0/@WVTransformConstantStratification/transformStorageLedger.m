function ledger = transformStorageLedger(self)
% Return known transform storage and explicitly opaque internal storage.
%
% This hidden benchmark contract excludes canonical model coefficients and
% forcing state. It reports transform-owned mappings, buffers, and vertical
% matrices exactly while keeping MATLAB-internal FFT storage opaque.
%
% - Topic: Inspect transform storage
% - Declaration: ledger = transformStorageLedger()
% - Returns ledger: exact application-owned arrays, opaque records, and aggregate byte counts
% - Developer: true

entries = self.fastTransform.storageLedger();
matrixNames = ["DCT" "iDCT" "DST" "iDST"];
for matrixName = matrixNames
    value = self.(matrixName);
    info = whos("value");
    entries(end+1,1) = ledgerEntry("vertical.matrix." + matrixName,"WVGeometryDoublyPeriodicStratifiedConstant","vertical-matrix","Dense vertical transform operator",string(class(value)),double(size(value)),double(info.bytes),"persistent","allocated","exact","real",double(info.bytes)); %#ok<AGROW>
end

isExact = string({entries.byteStatus}) == "exact";
isAllocated = string({entries.allocationState}) == "allocated";
isPersistent = string({entries.persistence}) == "persistent";
isTransient = string({entries.persistence}) == "transient";
knownPersistentBytes = sum([entries(isExact & isAllocated & isPersistent).bytes]);
knownTransientBytes = [entries(isExact & isAllocated & isTransient).bytes];
maximumKnownTransientBytes = max([0 knownTransientBytes]);
ledger = struct( ...
    "schema","transform-storage-ledger-v1", ...
    "scope","explicit-transform-storage", ...
    "implementation","builtin", ...
    "fourierStorageType",self.fastTransform.fourierStorageLayout.fourierStorageType, ...
    "entries",entries, ...
    "knownPersistentBytes",knownPersistentBytes, ...
    "knownTransientBytes",sum(knownTransientBytes), ...
    "maximumKnownTransientBytes",maximumKnownTransientBytes, ...
    "knownMaximumLiveBytes",knownPersistentBytes+maximumKnownTransientBytes, ...
    "opaqueEntryCount",nnz(string({entries.byteStatus}) == "opaque"), ...
    "hasPersistentFullSpectrum",any(isAllocated & isPersistent & string({entries.category}) == "spectrum-buffer" & string({entries.storageType}) == "full-complex"), ...
    "matlabInternalStorage","unresolved");
end

function value = ledgerEntry(identifier,owner,category,purpose,className,shape,bytes,persistence,allocationState,byteStatus,storageType,potentialBytes)
value = struct("identifier",identifier,"owner",owner,"category",category,"purpose",purpose,"className",className,"shape",double(shape),"bytes",double(bytes),"persistence",persistence,"allocationState",allocationState,"byteStatus",byteStatus,"storageType",storageType,"potentialBytes",double(potentialBytes));
end
