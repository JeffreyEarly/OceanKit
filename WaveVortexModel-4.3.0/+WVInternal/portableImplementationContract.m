function contract = portableImplementationContract(actualType,typeIdentifier,capabilityStatus,reason,payload)
% Build and validate a paired MATLAB/C++ implementation contract.
arguments
    actualType (1,1) string
    typeIdentifier (1,1) string
    capabilityStatus (1,1) string {mustBeMember(capabilityStatus,["supported","unavailable","versionMismatch","invalidContract"])}
    reason (1,1) string
    payload (1,1) struct
end

if strlength(actualType) == 0 || strlength(typeIdentifier) == 0
    capabilityStatus = "invalidContract";
    reason = "Portable type identifiers must be nonempty.";
    payload = struct();
elseif capabilityStatus == "supported" && actualType ~= typeIdentifier
    capabilityStatus = "invalidContract";
    reason = "An inherited portable contract cannot advertise a different MATLAB class.";
    payload = struct();
    typeIdentifier = actualType;
end

contract = struct( ...
    "schemaIdentifier","wave-vortex-portable-pair-v1", ...
    "schemaVersion",uint32(1), ...
    "typeIdentifier",typeIdentifier, ...
    "contractVersion",uint32(1), ...
    "capabilityStatus",capabilityStatus, ...
    "reason",reason, ...
    "payload",payload);
end
