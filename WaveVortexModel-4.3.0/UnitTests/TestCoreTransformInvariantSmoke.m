classdef TestCoreTransformInvariantSmoke < matlab.unittest.TestCase
    methods (Test, TestTags = "smoke")
        function representativeInvariantPath(testCase)
            wvt = WVTransformConstantStratification([6e3 4e3 1e3],[8 6 9],N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=true);
            seedRandomNumberGenerator(testCase,14001);
            wvt.initWithRandomFlow(uvMax=0.01);

            spectralIndices = (1:prod(wvt.spectralMatrixSize))';
            [kMode,lMode,jMode] = wvt.modeNumberFromIndex(spectralIndices);
            testCase.verifyEqual(wvt.indexFromModeNumber(kMode,lMode,jMode),spectralIndices)

            [u,v,w,eta] = wvt.variableWithName("u","v","w","eta");
            [Ap,Am,A0] = wvt.transformUVWEtaToWaveVortex(u,v,w,eta);
            testCase.verifyRelative(Ap,wvt.Ap,2e-12,"constant nonhydrostatic coefficient round trip")
            testCase.verifyRelative(Am,wvt.Am,2e-12,"constant nonhydrostatic coefficient round trip")
            testCase.verifyRelative(A0,wvt.A0,2e-12,"constant nonhydrostatic coefficient round trip")
            testCase.verifyRelative(wvt.totalEnergySpatiallyIntegrated,wvt.totalEnergy,2e-12,"constant nonhydrostatic Parseval energy")

            [X,~,Z] = wvt.xyzGrid;
            field = cos(4*pi*X/wvt.Lx).*cos(2*pi*(Z+wvt.Lz)/wvt.Lz);
            testCase.verifyRelative(wvt.diffX(field),-(4*pi/wvt.Lx)*sin(4*pi*X/wvt.Lx).*cos(2*pi*(Z+wvt.Lz)/wvt.Lz),2e-12,"constant nonhydrostatic horizontal derivative")
            testCase.verifyRelative(wvt.diffZF(field),-(2*pi/wvt.Lz)*cos(4*pi*X/wvt.Lx).*sin(2*pi*(Z+wvt.Lz)/wvt.Lz),2e-12,"constant nonhydrostatic vertical derivative")

            resized = wvt.waveVortexTransformWithResolution([9 7 8]);
            testCase.verifyEqual(resized.t,wvt.t)
            explicit = wvt.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyFalse(explicit.shouldAntialias)
            testCase.verifyTrue(explicit.hasForcingWithName("antialias filter"))

            for horizontalSize = {[8 7],[9 6]}
                geometry = WVGeometryDoublyPeriodic([6e3 4e3],horizontalSize{1},shouldAntialias=false);
                [X,Y] = ndgrid(geometry.x,geometry.y);
                parityField = cos(4*pi*X/geometry.Lx+2*pi*Y/geometry.Ly);
                parityCoefficients = geometry.transformFromSpatialDomainWithFourier(parityField);
                reconstructed = geometry.transformToSpatialDomainWithFourier(parityCoefficients);
                testCase.verifyRelative(reconstructed,parityField,2e-12,"mixed-parity Fourier round trip")
            end
        end
    end

    methods (Access = private)
        function verifyRelative(testCase,actual,expected,tolerance,diagnostic)
            relativeError = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
            testCase.verifyLessThanOrEqual(relativeError,tolerance,diagnostic)
        end
    end
end
