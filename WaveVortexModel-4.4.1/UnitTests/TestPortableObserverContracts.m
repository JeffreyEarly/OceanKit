classdef TestPortableObserverContracts < matlab.unittest.TestCase
    properties
        model
    end

    methods (TestMethodSetup)
        function createModel(testCase)
            wvt = WVTransformConstantStratification([1000 1000 100],[8 6 5],N0=5.2e-3,latitude=45,isHydrostatic=false);
            testCase.model = WVModel(wvt);
        end
    end

    methods (Test, TestTags={'full'})
        function baseClassIsUnavailable(testCase)
            observer = WVObservingSystem(testCase.model,"base observer");
            contract = observer.portableImplementationContract();
            testCase.verifyEqual(contract.capabilityStatus,"unavailable");
            testCase.verifyEqual(contract.typeIdentifier,"WVObservingSystem");
        end

        function supportedObserversReturnExactVersionedContracts(testCase)
            phi = zeros(testCase.model.wvt.Nx,testCase.model.wvt.Ny,testCase.model.wvt.Nz);
            observers = {
                WVCoefficients(testCase.model), ...
                WVEulerianFields(testCase.model,fieldNames="u"), ...
                WVMooring(testCase.model,x=0,y=0,trackedFieldNames="u"), ...
                WVLagrangianParticles(testCase.model,name="particles",x=0,y=0,z=0,trackedFieldNames="u"), ...
                WVTracer(testCase.model,name="tracer",phi=phi)};
            expected = ["WVCoefficients" "WVEulerianFields" "WVMooring" "WVLagrangianParticles" "WVTracer"];
            for iObserver = 1:numel(observers)
                contract = observers{iObserver}.portableImplementationContract();
                testCase.verifyEqual(contract.schemaIdentifier,"wave-vortex-portable-pair-v1");
                testCase.verifyEqual(contract.schemaVersion,uint32(1));
                testCase.verifyEqual(contract.typeIdentifier,expected(iObserver));
                testCase.verifyEqual(contract.contractVersion,uint32(1));
                testCase.verifyEqual(contract.capabilityStatus,"supported");
                testCase.verifyEqual(contract.reason,"");
                testCase.verifyTrue(isscalar(contract));
            end
        end

        function inheritedContractCannotAdvertiseParentClass(testCase)
            phi = zeros(testCase.model.wvt.Nx,testCase.model.wvt.Ny,testCase.model.wvt.Nz);
            observer = WVTestUnregisteredTracer(testCase.model,phi=phi);
            contract = observer.portableImplementationContract();
            testCase.verifyEqual(contract.capabilityStatus,"invalidContract");
            testCase.verifyEqual(contract.typeIdentifier,"WVTestUnregisteredTracer");
        end

        function testPairUsesDistinctIdentity(testCase)
            phi = zeros(testCase.model.wvt.Nx,testCase.model.wvt.Ny,testCase.model.wvt.Nz);
            observer = WVTestPortableTracer(testCase.model,phi=phi);
            contract = observer.portableImplementationContract();
            testCase.verifyEqual(contract.capabilityStatus,"supported");
            testCase.verifyEqual(contract.typeIdentifier,"WVTestPortableTracer");
            testCase.verifyEqual(contract.payload.stateShape,uint64(size(phi)));
        end

        function pointDiagnosticUsesDistinctConfigurationAndSampling(testCase)
            observer = WVTestPortablePointDiagnostic(testCase.model,x=[10; 20],y=[30; 40],z=[-10; -20],fieldName="u",scale=2.5,offset=-1.25,interpolation="spline");
            contract = observer.portableImplementationContract();
            testCase.verifyEqual(contract.capabilityStatus,"supported");
            testCase.verifyEqual(contract.typeIdentifier,"WVTestPortablePointDiagnostic");
            testCase.verifyEqual(contract.payload.fieldNames,"u");
            testCase.verifyEqual(contract.payload.x,[10; 20]);
            testCase.verifyEqual(contract.payload.outputScale,2.5);
            testCase.verifyEqual(contract.payload.outputOffset,-1.25);
            testCase.verifyEqual(contract.payload.trackedFieldInterpolation,"spline");
        end
    end
end
