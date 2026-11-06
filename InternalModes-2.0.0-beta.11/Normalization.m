classdef Normalization
    % Enumerate legacy v1 internal-mode scaling conventions.
    %
    % This enumeration is retained for the v1 implementation. V2 basis sets
    % accept text names such as "rmsF" and "waveEnergy" and reject these enum
    % values. The v2 catalogue is reported by normalizationNames on each basis.
    %
    % - Topic: Configure normalization
    % - Declaration: classdef Normalization
    enumeration
        unity, kConstant, uMax, wMax, surfacePressure, geostrophic, depth, generalizedPotentialEnstrophy
    end
end
