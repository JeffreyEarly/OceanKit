classdef IMInternalModesBasisTestAccess < IMInternalModesBasis
    methods
        function self = IMInternalModesBasisTestAccess(source)
            arguments
                source (1,1) IMInternalModesBasis
            end
            self@IMInternalModesBasis(solver=source.solver,evp=source.evp,nativeModes=source.nativeModes, ...
                eigenvalues=source.eigenvalues,modeNumber=source.modeNumber, ...
                modeSelectionDiagnostics=source.modeSelectionDiagnostics,normalization=source.normalization,metadata=source.metadata);
        end

        function mask = classifyTransformColumns(self,sampled,targetGram,targetMajorantGram)
            mask = self.internalModesTransformActiveMask("F",sampled,targetGram,targetMajorantGram);
        end
    end
end
