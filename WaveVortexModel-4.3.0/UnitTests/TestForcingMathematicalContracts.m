classdef TestForcingMathematicalContracts < matlab.unittest.TestCase
    properties
        repositoryRoot (1,1) string
    end

    methods (TestMethodSetup)
        function locateRepository(testCase)
            testCase.repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test, TestTags="full")
        function pseudoTopographicDocumentationDescribesProjection(testCase)
            page = testCase.generatedPage("forcing/wvpseudotopographicwavegeneration/index.md");
            requiredText = [
                "\mathbf{U}_{\mathrm{bt}}"
                "\widehat{g}_b"
                "Apm_TE_factor"
                "bottom pressure work"
                "conjPhase"
                "`phase` for `Am`"
                "incoming `A0` tendency is unchanged"
                "WVAdaptiveDamping"
                ];
            for token = requiredText'
                testCase.verifySubstring(page,token);
            end
            testCase.verifyFalse(contains(page,"\boldsymbol"));
        end

        function betaPlaneAnalyticalSignIsNegativeBetaV(testCase)
            wvt = WVTransformBarotropicQG([40e3 30e3],[8 6],latitude=45,shouldAntialias=false);
            wvt.removeAllForcing();
            amplitude = 0.02;
            k = 2*pi/wvt.Lx;
            streamfunction = @(x,y,z) amplitude*cos(k*x);
            wvt.initWithGeostrophicStreamfunction(streamfunction);

            forcing = WVBetaPlanePVAdvection(wvt);
            tendency = forcing.addPotentialVorticitySpatialForcing(wvt,zeros(wvt.spatialMatrixSize));
            expectedV = -amplitude*k*sin(k*wvt.X);
            expectedTendency = -wvt.beta*expectedV;
            testCase.verifyEqual(tendency,expectedTendency,AbsTol=1e-14*max(1,max(abs(expectedTendency),[],"all")));
        end

        function betaPlanePhysicalPathMatchesQGPVContract(testCase)
            transforms = TestForcingMathematicalContracts.qgTransforms();
            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                wvt.removeAllForcing();
                wvt.A0 = TestForcingMathematicalContracts.deterministicCoefficients(size(wvt.A0)).*wvt.geostrophicComponent.maskA0;
                forcing = WVBetaPlanePVAdvection(wvt);
                tendency = forcing.addPotentialVorticitySpatialForcing(wvt,zeros(wvt.spatialMatrixSize));
                testCase.verifyEqual(tendency,-wvt.beta*wvt.v,AbsTol=1e-14*max(1,max(abs(wvt.beta*wvt.v),[],"all")));
            end
        end

        function betaPlaneSpectralPathAffectsOnlyGeostrophicA0(testCase)
            transforms = TestForcingMathematicalContracts.waveTransforms();
            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                wvt.removeAllForcing();
                coefficients = TestForcingMathematicalContracts.deterministicCoefficients(size(wvt.A0));
                wvt.A0 = coefficients.*wvt.geostrophicComponent.maskA0+0.25*coefficients.*wvt.mdaComponent.maskA0;
                wvt.Ap = TestForcingMathematicalContracts.deterministicCoefficients(size(wvt.Ap)).*wvt.waveComponent.maskAp;
                wvt.Am = conj(TestForcingMathematicalContracts.deterministicCoefficients(size(wvt.Am))).*wvt.waveComponent.maskAm;

                forcing = WVBetaPlanePVAdvection(wvt);
                incomingFp = 2*ones(size(wvt.Ap));
                incomingFm = -3*ones(size(wvt.Am));
                incomingF0 = 5*ones(size(wvt.A0));
                [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,incomingFp,incomingFm,incomingF0);
                contribution = F0-incomingF0;

                testCase.verifyEqual(Fp,incomingFp);
                testCase.verifyEqual(Fm,incomingFm);
                testCase.verifyEqual(contribution(~logical(wvt.geostrophicComponent.maskA0)),zeros(nnz(~logical(wvt.geostrophicComponent.maskA0)),1));
                testCase.verifyEqual(contribution(logical(wvt.mdaComponent.maskA0)),zeros(nnz(wvt.mdaComponent.maskA0),1));
                testCase.verifyEqual(contribution(wvt.Kh == 0),zeros(nnz(wvt.Kh == 0),1));

                qgpvTendency = wvt.transformToSpatialDomainWithF(A0=wvt.A0_QGPV_factor.*contribution);
                geostrophicV = wvt.transformToSpatialDomainWithF(A0=wvt.VA0.*wvt.A0);
                scale = max(1,max(abs(wvt.beta*geostrophicV),[],"all"));
                testCase.verifyEqual(qgpvTendency,-wvt.beta*geostrophicV,AbsTol=1e-13*scale);
            end
        end

        function stratifiedQGVerticalDiffusivityMatchesInducedQGPVTendency(testCase)
            wvt = TestForcingMathematicalContracts.stratifiedQGTransform();
            wvt.removeAllForcing();
            wvt.A0 = TestForcingMathematicalContracts.deterministicCoefficients(size(wvt.A0)).*wvt.geostrophicComponent.maskA0;
            withCorrection = WVVerticalDiffusivity(wvt,kappa_z=7e-6,shouldForceMeanDensityAnomaly=true);
            withoutCorrection = WVVerticalDiffusivity(wvt,kappa_z=7e-6,shouldForceMeanDensityAnomaly=false);

            actual = withCorrection.addPotentialVorticitySpatialForcing(wvt,zeros(wvt.spatialMatrixSize));
            displacementTendency = withCorrection.kappa_z*wvt.diffZG(wvt.eta,n=2);
            expected = -wvt.f*wvt.diffZG(displacementTendency,n=1);
            scale = max(1,max(abs(expected),[],"all"));
            testCase.verifyEqual(actual,expected,AbsTol=1e-13*scale);
            testCase.verifyEqual(withoutCorrection.addPotentialVorticitySpatialForcing(wvt,zeros(wvt.spatialMatrixSize)),actual);
            testCase.verifyFalse(isa(wvt,"WVMeanDensityAnomalyMethods"));
            testCase.verifyEqual(string({wvt.primaryFlowComponents.shortName}),"geostrophic");
            testCase.verifyEqual(nnz(wvt.geostrophicComponent.maskA0(wvt.Kh == 0)),0);
        end

        function meanDensityAnomalyFlagControlsOnlyUniformMDASource(testCase)
            Lxyz = [40e3 30e3 2e3];
            Nxyz = [8 6 5];
            N2 = @(z)2e-5*exp(z/4000);
            transforms = {
                WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false)
                WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false)
                };
            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                withCorrection = WVVerticalDiffusivity(wvt,kappa_z=7e-6,shouldForceMeanDensityAnomaly=true);
                withoutCorrection = WVVerticalDiffusivity(wvt,kappa_z=7e-6,shouldForceMeanDensityAnomaly=false);
                zeroField = zeros(wvt.spatialMatrixSize);
                if wvt.isHydrostatic
                    [~,~,withTendency] = withCorrection.addHydrostaticSpatialForcing(wvt,zeroField,zeroField,zeroField);
                    [~,~,withoutTendency] = withoutCorrection.addHydrostaticSpatialForcing(wvt,zeroField,zeroField,zeroField);
                else
                    [~,~,~,withTendency] = withCorrection.addNonhydrostaticSpatialForcing(wvt,zeroField,zeroField,zeroField,zeroField);
                    [~,~,~,withoutTendency] = withoutCorrection.addNonhydrostaticSpatialForcing(wvt,zeroField,zeroField,zeroField,zeroField);
                end
                expectedCorrection = -withCorrection.kappa_z*(zeroField+shiftdim(wvt.dLnN2,-2));
                testCase.verifyEqual(withCorrection.dLnN2,shiftdim(wvt.dLnN2,-2));
                testCase.verifyEqual(withTendency-withoutTendency,expectedCorrection,AbsTol=10*eps);

                if wvt.isHydrostatic
                    [Fp,Fm,F0] = wvt.transformUVEtaToWaveVortex(zeroField,zeroField,expectedCorrection);
                else
                    [Fp,Fm,F0] = wvt.transformUVWEtaToWaveVortex(zeroField,zeroField,zeroField,expectedCorrection);
                end
                scale = max(1,max(abs(F0),[],"all"));
                testCase.verifyEqual(Fp,zeros(size(Fp)),AbsTol=1e-14*scale);
                testCase.verifyEqual(Fm,zeros(size(Fm)),AbsTol=1e-14*scale);
                testCase.verifyEqual(F0(~logical(wvt.mdaComponent.maskA0)),zeros(nnz(~logical(wvt.mdaComponent.maskA0)),1),AbsTol=1e-14*scale);
                testCase.verifyGreaterThan(max(abs(F0(logical(wvt.mdaComponent.maskA0)))),0);
            end
        end

        function correctionOptionIsNaturallyInertWhenUnavailable(testCase)
            constant = WVTransformConstantStratification([40e3 30e3 2e3],[8 6 5],N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false);
            withCorrection = WVVerticalDiffusivity(constant,shouldForceMeanDensityAnomaly=true);
            withoutCorrection = WVVerticalDiffusivity(constant,shouldForceMeanDensityAnomaly=false);
            zeroField = zeros(constant.spatialMatrixSize);
            [~,~,withTendency] = withCorrection.addHydrostaticSpatialForcing(constant,zeroField,zeroField,zeroField);
            [~,~,withoutTendency] = withoutCorrection.addHydrostaticSpatialForcing(constant,zeroField,zeroField,zeroField);
            testCase.verifyEqual(withCorrection.dLnN2,0);
            testCase.verifyEqual(withTendency,withoutTendency);

            barotropic = WVTransformBarotropicQG([40e3 30e3],[8 6],latitude=45,shouldAntialias=false);
            forcing = testCase.verifyWarningFree(@()WVVerticalDiffusivity(barotropic,shouldForceMeanDensityAnomaly=true));
            testCase.verifyEqual(forcing.dLnN2,0);
            testCase.verifyError(@()barotropic.addForcing(forcing),'');
        end

        function generatedDocumentationStatesVerifiedContracts(testCase)
            beta = testCase.generatedPage("forcing/wvbetaplanepvadvection/index.md");
            vertical = testCase.generatedPage("forcing/closures/wvverticaldiffusivity/index.md");
            guide = string(fileread(fullfile(testCase.repositoryRoot,"docs","users-guide","adding-forcing.md")));
            testCase.verifySubstring(beta,"-\beta v_g");
            testCase.verifySubstring(beta,"not a full");
            testCase.verifySubstring(beta,"beta-plane treatment of internal-wave dynamics");
            testCase.verifyFalse(contains(beta,"I suspect"));
            testCase.verifySubstring(vertical,"-f\kappa_z\partial_{zzz}\eta");
            testCase.verifySubstring(vertical,"not a mean-density-anomaly component");
            testCase.verifySubstring(vertical,"does not modify the wave modes");
            testCase.verifySubstring(guide,"conservation of $$q+\beta y$$");
        end
    end

    methods (Access=private)
        function page = generatedPage(testCase,relativePath)
            page = string(fileread(fullfile(testCase.repositoryRoot,"docs","classes",relativePath)));
        end
    end

    methods (Static, Access=private)
        function transforms = waveTransforms()
            Lxyz = [40e3 30e3 2e3];
            Nxyz = [8 6 5];
            N2 = @(z)2e-5*exp(z/4000);
            transforms = {
                WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false)
                WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false)
                WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false)
                WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false)
                };
        end

        function transforms = qgTransforms()
            transforms = {
                TestForcingMathematicalContracts.stratifiedQGTransform()
                WVTransformBarotropicQG([40e3 30e3],[8 6],latitude=45,shouldAntialias=false)
                };
        end

        function wvt = stratifiedQGTransform()
            wvt = WVTransformStratifiedQG([40e3 30e3 2e3],[8 6 5],N2=@(z)2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
        end

        function coefficients = deterministicCoefficients(matrixSize)
            values = reshape(1:prod(matrixSize),matrixSize);
            coefficients = complex(values,reshape(prod(matrixSize):-1:1,matrixSize));
        end
    end
end
