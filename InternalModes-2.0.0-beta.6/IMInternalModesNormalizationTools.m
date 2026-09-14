classdef (Hidden, Sealed) IMInternalModesNormalizationTools
    % Define normalization availability and defaults for internal-mode EVPs.

    methods (Static)
        function names = available(evp)
            names = ["unity","uMax","wMax"];
            if evp.modeFamily ~= "meanDensityAnomaly"
                names(end+1) = "surfacePressure";
            end
            if evp.modeFamily == "hydrostatic"
                names(end+1) = "geostrophic";
                if evp.innerProduct("F").hasInnerProduct
                    names(end+1) = "depth";
                end
            end
            if string(evp.name) == "waveModesAtWavenumber"
                names(end+1) = "kConstant";
            end
            if string(evp.name) == "geostrophicGeneralizedPotentialEnstrophyModes"
                names(end+1) = "generalizedPotentialEnstrophy";
            end
        end

        function name = default(evp)
            if string(evp.name) == "geostrophicGeneralizedPotentialEnstrophyModes"
                name = "generalizedPotentialEnstrophy";
            elseif string(evp.name) == "geostrophicAPVModes"
                name = "depth";
            elseif evp.modeFamily == "hydrostatic"
                name = "geostrophic";
            elseif string(evp.name) == "waveModesAtWavenumber"
                name = "kConstant";
            else
                name = "unity";
            end
        end
    end
end
