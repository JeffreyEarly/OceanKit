classdef IMOrientationRawGOverrideBasis < IMInternalModesBasis
    methods
        function self = IMOrientationRawGOverrideBasis(source)
            arguments
                source (1,1) IMInternalModesBasis
            end
            self@IMInternalModesBasis(solver=source.solver,evp=source.evp,nativeModes=source.nativeModes, ...
                eigenvalues=source.eigenvalues,modeNumber=source.modeNumber, ...
                modeSelectionDiagnostics=source.modeSelectionDiagnostics,normalization=source.normalization,metadata=source.metadata);
        end

        function values = rawVariable(self,variable,z)
            values = rawVariable@IMInternalModesBasis(self,variable,z);
            if string(variable) == "G"
                values = -values;
            end
        end
    end
end
