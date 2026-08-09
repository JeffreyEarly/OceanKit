classdef TestInertialOscillationMethods < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[15e3, 15e3, 1300]);
        Nxyz = struct('Nx8Ny8Nz30',[8 8 30]);
        transform = {'constant-hydrostatic','constant-nonhydrostatic','hydrostatic','boussinesq'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant-hydrostatic'
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,isHydrostatic=true,shouldAntialias=false);
                case 'constant-nonhydrostatic'
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,isHydrostatic=false,shouldAntialias=false);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=@(z) (5.2e-3)^2*ones(size(z)),shouldAntialias=false);
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=@(z) (5.2e-3)^2*ones(size(z)),shouldAntialias=false);
            end
            % testCase.wvt.addOperation(testCase.wvt.operationForDynamicalVariable('u','v','eta','w',flowComponent=testCase.wvt.flowComponent('inertial')));
            testCase.solutionGroup = WVInertialOscillationComponent(testCase.wvt);
        end
    end

    methods (TestMethodSetup)
        function resetTransform(self)
            self.wvt.removeAll();
            self.wvt.t = 0;
        end
    end

    methods (Test, TestTags = "full")
        function testRemoveAllInertialMotions(self)
            seedRandomNumberGenerator(self,21517);
            % In this test we intialize with a random flow state, confirm
            % that both total energy and inertial energy are present,
            % remove all the inertial energy, and then confirm that there
            % is no inertial energy remaining, and that the total energy is
            % the same as the initial total, minus the initial inertial.
            self.wvt.initWithRandomFlow(uvMax=0.02);

            initialTotalEnergy = self.wvt.totalEnergy;
            initialInertialEnergy = self.wvt.inertialEnergy;
            self.verifyGreaterThan(initialTotalEnergy,0.0);
            self.verifyGreaterThan(initialInertialEnergy,0.0);

            self.wvt.removeAllInertialMotions();
            finalTotalEnergy = self.wvt.totalEnergy;
            finalInertialEnergy = self.wvt.inertialEnergy;

            self.verifyEqual(finalInertialEnergy,0.0);
            self.verifyEqual(finalTotalEnergy,initialTotalEnergy-initialInertialEnergy, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testInitWithInertialMotions(self)
            seedRandomNumberGenerator(self,37547);
            % In this test we expect all existing motions to be removed
            % when we initialize.
            %
            % VERY IMPORTANT: We do *not* get back exactly what we put in
            % because we are de-aliasing the signal. Hence, this experiment
            % was designed assuming Nz=30 (and thus Nj=20).
            U_io = 0.2;
            Ld = self.wvt.Lz/2;
            theta = pi/3;
            u_NIO = @(z) U_io*cos(theta)*exp((z/Ld));
            v_NIO = @(z) U_io*sin(theta)*exp((z/Ld));

            % Populate the flow field with junk...
            self.wvt.initWithRandomFlow(uvMax=0.02);

            % call our initWithInertialMotions method
            self.wvt.initWithInertialMotions(u_NIO,v_NIO);

            % now verify that only inertial oscillations are part of the
            % solution.
            self.verifyThat(u_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.u,relTol=1e-3),'u_tot');
            self.verifyThat(v_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.v,relTol=1e-3),'v_tot');
        end

        function testSetInertialMotions(self)
            seedRandomNumberGenerator(self,58451);
            U_io = 0.2;
            Ld = self.wvt.Lz/2;
            theta = 0;
            u_NIO = @(z) U_io*cos(theta)*exp((z/Ld));
            v_NIO = @(z) U_io*sin(theta)*exp((z/Ld));

            % Populate the flow field with junk...
            self.wvt.initWithRandomFlow(uvMax=0.02);

            initialTotalEnergy = self.wvt.totalEnergy;
            initialInertialEnergy = self.wvt.inertialEnergy;

            % Now overwrite *only* the inertial stuff, other stuff should
            % remain.
            self.wvt.setInertialMotions(u_NIO,v_NIO);

            finalTotalEnergy = self.wvt.totalEnergy;
            finalInertialEnergy = self.wvt.inertialEnergy;

            % Now confirm that the inertial solution matches AND that the
            % original non-inertial energy stayed the same.
            self.verifyThat(u_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.u_io,relTol=1e-3),'u_io');
            self.verifyThat(v_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.v_io,relTol=1e-3),'v_io');
            self.verifyEqual(finalTotalEnergy-finalInertialEnergy,initialTotalEnergy-initialInertialEnergy, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testRandomFlowProducesConjugateInertialCoefficients(self)
            seedRandomNumberGenerator(self,74653);
            self.wvt.initWithRandomFlow(uvMax=0.02);

            inertialApMask = logical(self.wvt.inertialComponent.maskAp);
            inertialAmMask = logical(self.wvt.inertialComponent.maskAm);
            self.verifyGreaterThan(nnz(self.wvt.Ap(inertialApMask)),0)
            self.verifyGreaterThan(nnz(self.wvt.Am(inertialAmMask)),0)
            self.verifyEqual(self.wvt.Am(inertialAmMask),conj(self.wvt.Ap(inertialApMask)))

            primaryWaveApMask = logical(self.wvt.waveComponent.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap));
            primaryWaveAmMask = logical(self.wvt.waveComponent.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Am));
            self.verifyGreaterThan(nnz(self.wvt.Ap(primaryWaveApMask)),0)
            self.verifyGreaterThan(nnz(self.wvt.Am(primaryWaveAmMask)),0)
        end

        function testAddInertialMotionsSuperposesValidState(self)
            seedRandomNumberGenerator(self,86477);
            [u1,v1,u2,v2] = TestInertialOscillationMethods.inertialProfiles(self.wvt);
            self.wvt.initWithRandomFlow(uvMax=0.02);

            self.wvt.setInertialMotions(u2,v2);
            inertialApMask = logical(self.wvt.inertialComponent.maskAp);
            inertialAmMask = logical(self.wvt.inertialComponent.maskAm);
            incrementAp = zeros(size(self.wvt.Ap));
            incrementAm = zeros(size(self.wvt.Am));
            incrementAp(inertialApMask) = self.wvt.Ap(inertialApMask);
            incrementAm(inertialAmMask) = self.wvt.Am(inertialAmMask);
            incrementU = self.wvt.u_io;
            incrementV = self.wvt.v_io;

            self.wvt.setInertialMotions(u1,v1);
            beforeAp = self.wvt.Ap;
            beforeAm = self.wvt.Am;
            beforeA0 = self.wvt.A0;
            beforeUInertial = self.wvt.u_io;
            beforeVInertial = self.wvt.v_io;
            beforeU = self.wvt.u;
            beforeV = self.wvt.v;

            expectedAp = beforeAp + incrementAp;
            expectedAm = beforeAm + incrementAm;
            expectedTotalEnergy = TestInertialOscillationMethods.totalEnergyForCoefficients(self.wvt,expectedAp,expectedAm,beforeA0);
            expectedInertialEnergy = TestInertialOscillationMethods.inertialEnergyForCoefficients(self.wvt,expectedAp,expectedAm);

            self.wvt.addInertialMotions(u2,v2);

            self.verifyEqual(self.wvt.Ap,expectedAp,AbsTol=100*eps)
            self.verifyEqual(self.wvt.Am,expectedAm,AbsTol=100*eps)
            self.verifyEqual(self.wvt.A0,beforeA0)
            self.verifyEqual(self.wvt.Ap(~inertialApMask),beforeAp(~inertialApMask))
            self.verifyEqual(self.wvt.Am(~inertialAmMask),beforeAm(~inertialAmMask))
            self.verifyEqual(self.wvt.u_io,beforeUInertial+incrementU,AbsTol=100*eps)
            self.verifyEqual(self.wvt.v_io,beforeVInertial+incrementV,AbsTol=100*eps)
            self.verifyEqual(self.wvt.u,beforeU+incrementU,AbsTol=100*eps)
            self.verifyEqual(self.wvt.v,beforeV+incrementV,AbsTol=100*eps)
            self.verifyEqual(self.wvt.totalEnergy,expectedTotalEnergy,AbsTol=100*eps,RelTol=100*eps)
            self.verifyEqual(self.wvt.inertialEnergy,expectedInertialEnergy,AbsTol=100*eps,RelTol=100*eps)
        end

        function testAddInertialMotionsAfterRandomFlow(self)
            seedRandomNumberGenerator(self,104729);
            [u1,v1] = TestInertialOscillationMethods.inertialProfiles(self.wvt);

            self.wvt.initWithInertialMotions(u1,v1);
            Ap1 = self.wvt.Ap;
            Am1 = self.wvt.Am;
            A01 = self.wvt.A0;
            uInertial1 = self.wvt.u_io;
            vInertial1 = self.wvt.v_io;
            u1Total = self.wvt.u;
            v1Total = self.wvt.v;

            self.wvt.initWithRandomFlow(uvMax=0.02);
            Ap2 = self.wvt.Ap;
            Am2 = self.wvt.Am;
            A02 = self.wvt.A0;
            uInertial2 = self.wvt.u_io;
            vInertial2 = self.wvt.v_io;
            u2Total = self.wvt.u;
            v2Total = self.wvt.v;
            inertialApMask = logical(self.wvt.inertialComponent.maskAp);
            inertialAmMask = logical(self.wvt.inertialComponent.maskAm);
            self.verifyEqual(Am2(inertialAmMask),conj(Ap2(inertialApMask)))

            expectedAp = Ap1 + Ap2;
            expectedAm = Am1 + Am2;
            expectedA0 = A01 + A02;
            expectedTotalEnergy = TestInertialOscillationMethods.totalEnergyForCoefficients(self.wvt,expectedAp,expectedAm,expectedA0);
            expectedInertialEnergy = TestInertialOscillationMethods.inertialEnergyForCoefficients(self.wvt,expectedAp,expectedAm);

            self.wvt.addInertialMotions(u1,v1);

            self.verifyEqual(self.wvt.Ap,expectedAp,AbsTol=100*eps)
            self.verifyEqual(self.wvt.Am,expectedAm,AbsTol=100*eps)
            self.verifyEqual(self.wvt.A0,expectedA0)
            self.verifyEqual(self.wvt.Ap(~inertialApMask),Ap2(~inertialApMask))
            self.verifyEqual(self.wvt.Am(~inertialAmMask),Am2(~inertialAmMask))
            self.verifyEqual(self.wvt.u_io,uInertial1+uInertial2,AbsTol=100*eps)
            self.verifyEqual(self.wvt.v_io,vInertial1+vInertial2,AbsTol=100*eps)
            self.verifyEqual(self.wvt.u,u1Total+u2Total,AbsTol=100*eps)
            self.verifyEqual(self.wvt.v,v1Total+v2Total,AbsTol=100*eps)
            self.verifyEqual(self.wvt.totalEnergy,expectedTotalEnergy,AbsTol=100*eps,RelTol=100*eps)
            self.verifyEqual(self.wvt.inertialEnergy,expectedInertialEnergy,AbsTol=100*eps,RelTol=100*eps)
        end

    end

    methods (Static, Access=private)
        function [u1,v1,u2,v2] = inertialProfiles(wvt)
            u1 = @(z) 0.12*exp(z/(wvt.Lz/2));
            v1 = @(z) 0.04*exp(z/(wvt.Lz/3));
            u2 = @(z) -0.03*cos(pi*(z+wvt.Lz)/wvt.Lz);
            v2 = @(z) 0.07*cos(2*pi*(z+wvt.Lz)/wvt.Lz);
        end

        function energy = totalEnergyForCoefficients(wvt,Ap,Am,A0)
            energy = sum(wvt.Apm_TE_factor(:).*(abs(Ap(:)).^2+abs(Am(:)).^2)+wvt.A0_TE_factor(:).*abs(A0(:)).^2);
        end

        function energy = inertialEnergyForCoefficients(wvt,Ap,Am)
            maskAp = wvt.inertialComponent.maskAp;
            maskAm = wvt.inertialComponent.maskAm;
            energy = sum(wvt.Apm_TE_factor(:).*(maskAp(:).*abs(Ap(:)).^2+maskAm(:).*abs(Am(:)).^2));
        end
    end

end
