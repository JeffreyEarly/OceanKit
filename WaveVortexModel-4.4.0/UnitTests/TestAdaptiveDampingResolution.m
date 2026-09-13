classdef TestAdaptiveDampingResolution < matlab.unittest.TestCase
    methods (Test,TestTags="full")
        function collapsedAndOrdinaryFiltersHaveDefinedLimits(testCase)
            wvt = TestAdaptiveDampingResolution.transform();
            damping = WVAdaptiveDamping(wvt);
            dkl = min(wvt.dk,wvt.dl);
            for horizontalMaximum = [0 dkl 3*dkl]
                for verticalMaximum = [0 1 3]
                    [Qkl,Qj,kCutoff,~,jCutoff] = damping.spectralVanishingViscosityFilter(horizontalMaximum,verticalMaximum);
                    testCase.verifyTrue(all(isfinite(Qkl),"all"));
                    testCase.verifyTrue(all(isfinite(Qj),"all"));
                    if horizontalMaximum <= dkl
                        testCase.verifyEqual(Qkl,double(wvt.Kh>=horizontalMaximum & wvt.Kh>0));
                    else
                        expected = exp(-((wvt.Kh-horizontalMaximum)./(wvt.Kh-kCutoff)).^2);
                        expected(wvt.Kh<kCutoff) = 0; expected(wvt.Kh>horizontalMaximum) = 1;
                        testCase.verifyEqual(Qkl,expected);
                    end
                    if verticalMaximum <= 1
                        testCase.verifyEqual(Qj,double(wvt.J>=verticalMaximum & wvt.J>0));
                    else
                        expected = exp(-((wvt.J-verticalMaximum)./(wvt.J-jCutoff)).^2);
                        expected(wvt.J<jCutoff) = 0; expected(wvt.J>verticalMaximum) = 1;
                        testCase.verifyEqual(Qj,expected);
                    end
                end
            end
        end

        function explicitHorizontalFilterDefinesSupportedFullRate(testCase)
            wvt = WVTransformConstantStratification([17000 17000 1000],[5 5 9],N0=5.2e-3,isHydrostatic=true,shouldAntialias=false);
            wvt.setForcing(WVAntialiasing(wvt,Nj=2));
            testCase.verifyEqual(pi/wvt.effectiveHorizontalGridResolution,wvt.dk,RelTol=1e-14);
            damping = WVAdaptiveDamping(wvt);
            testCase.verifyTrue(all(isfinite(damping.damp),"all"));
            testCase.verifyEqual(damping.damp(wvt.Kh==0 & wvt.J==0),0);
            unresolved = WVTransformConstantStratification([17000 17000 1000],[3 3 9],N0=5.2e-3,isHydrostatic=true,shouldAntialias=false);
            unresolved.setForcing(WVAntialiasing(unresolved,Nj=2));
            testCase.verifyEqual(unresolved.effectiveHorizontalGridResolution,Inf);
            testCase.verifyError(@()WVAdaptiveDamping(unresolved),"WVAdaptiveDamping:UnresolvedHorizontalWavenumbers");
        end

        function changingOnlyVerticalResolutionRefreshesDamping(testCase)
            wvt = TestAdaptiveDampingResolution.transform();
            damping = WVAdaptiveDamping(wvt);
            wvt.setForcing(damping);
            dx = [];
            for retainedModes = [3 1 2 4]
                wvt.setForcing([WVAntialiasing(wvt,Nj=retainedModes),damping]);
                testCase.verifyEqual(wvt.effectiveJMax,retainedModes-1);
                if ~isempty(dx), testCase.verifyEqual(wvt.effectiveHorizontalGridResolution,dx); end
                dx = wvt.effectiveHorizontalGridResolution;
                testCase.verifyTrue(all(isfinite(damping.damp),"all"));
                reference = WVAdaptiveDamping(wvt);
                testCase.verifyEqual(damping.damp,reference.damp);
                testCase.verifyEqual(damping.j_no_damp,reference.j_no_damp);
                testCase.verifyEqual(damping.j_damp,reference.j_damp);
                if retainedModes >= 3
                    [Qkl,Qj] = reference.spectralVanishingViscosityFilter(pi/dx,retainedModes-1);
                    xy = dx/(pi^2);
                    zz = (pi*pi*wvt.Lr2(retainedModes)/(dx^2))*xy;
                    expected = -xy*Qkl.*(wvt.K.^2+wvt.L.^2)-zz*Qj.*(1./wvt.Lr2);
                    testCase.verifyEqual(damping.damp,expected);
                end
            end
        end

        function zeroVerticalResolutionPreservesUniformLargeScaleFlow(testCase)
            wvt = TestAdaptiveDampingResolution.transform();
            wvt.A0 = zeros(wvt.spectralMatrixSize);
            index = find(wvt.J==0 & wvt.Kh>0 & wvt.Kh<1.01*min(wvt.dk,wvt.dl),1);
            testCase.assertNotEmpty(index);
            wvt.A0(index) = 1e-6;
            damping = WVAdaptiveDamping(wvt);
            wvt.setForcing([WVAntialiasing(wvt,Nj=1),damping]);
            testCase.verifyGreaterThan(wvt.uvMax,0);
            zero = zeros(wvt.spectralMatrixSize);
            [Fp,Fm,F0] = damping.addSpectralForcing(wvt,zero,zero,zero);
            testCase.verifyEqual(Fp,zero); testCase.verifyEqual(Fm,zero); testCase.verifyEqual(F0,zero);
            testCase.verifyEqual(damping.damp(wvt.Kh==0 & wvt.J==0),0);
        end
    end
    methods (Static,Access=private)
        function wvt = transform()
            wvt = WVTransformConstantStratification([17000 17000 1000],[8 8 9],N0=5.2e-3,isHydrostatic=true,shouldAntialias=false);
        end
    end
end
