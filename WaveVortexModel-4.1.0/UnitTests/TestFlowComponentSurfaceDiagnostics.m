classdef TestFlowComponentSurfaceDiagnostics < matlab.unittest.TestCase

    properties
        wvt
        balancedComponent
    end

    properties (ClassSetupParameter)
        transform = {'constant','hydrostatic','boussinesq'}
    end

    methods (TestClassSetup)
        function classSetup(testCase,transform)
            Lxyz = [1000 500 500];
            Nxyz = [8 8 5];
            N2 = @(z) (5.2e-3)^2*ones(size(z));

            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,latitude=33,isHydrostatic=false,shouldAntialias=false);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,shouldAntialias=false);
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,shouldAntialias=false);
            end

            rng(1);
            testCase.wvt.initWithRandomFlow(uvMax=0.01);
            testCase.wvt.t = 3600;

            testCase.balancedComponent = testCase.wvt.geostrophicComponent + testCase.wvt.mdaComponent;
            testCase.balancedComponent.name = 'balanced';
            testCase.balancedComponent.shortName = 'balanced';
            testCase.balancedComponent.abbreviatedName = 'balanced';
            testCase.wvt.addFlowComponent(testCase.balancedComponent);

            components = [testCase.wvt.geostrophicComponent,testCase.wvt.waveComponent,testCase.wvt.inertialComponent,testCase.wvt.mdaComponent,testCase.balancedComponent];
            for iComponent = 1:length(components)
                testCase.wvt.addOperation(testCase.wvt.operationForKnownVariable('pi',flowComponent=components(iComponent)));
            end
        end
    end

    methods (Test)
        function testRegistrationAddsStandardVariables(testCase)
            testCase.verifyTrue(ismember("balanced",testCase.wvt.flowComponentNames));
            expectedNames = ["u_balanced","v_balanced","w_balanced","eta_balanced","p_balanced","ssh_balanced","ssu_balanced","ssv_balanced"];
            testCase.verifyTrue(all(testCase.wvt.hasVariableWithName(expectedNames{:})));
        end

        function testBalancedSurfaceClosure(testCase)
            variables = {'ssh','ssu','ssv'};
            for iVariable = 1:length(variables)
                name = variables{iVariable};
                balanced = testCase.wvt.variableWithName([name '_balanced']);
                geostrophic = testCase.wvt.variableWithName([name '_g']);
                mda = testCase.wvt.variableWithName([name '_mda']);
                testCase.verifyEqual(balanced,geostrophic+mda,AbsTol=1e-12);
            end
        end

        function testTotalSurfaceClosure(testCase)
            variables = {'ssh','ssu','ssv'};
            for iVariable = 1:length(variables)
                name = variables{iVariable};
                total = testCase.wvt.variableWithName(name);
                geostrophic = testCase.wvt.variableWithName([name '_g']);
                mda = testCase.wvt.variableWithName([name '_mda']);
                wave = testCase.wvt.variableWithName([name '_w']);
                inertial = testCase.wvt.variableWithName([name '_io']);
                testCase.verifyEqual(total,geostrophic+mda+wave+inertial,AbsTol=1e-12);
            end
        end

        function testComponentPressureUsesComponentPressureHeight(testCase)
            suffixes = {'g','mda','w','io','balanced'};
            for iSuffix = 1:length(suffixes)
                suffix = suffixes{iSuffix};
                pressure = testCase.wvt.variableWithName(['p_' suffix]);
                pressureHeight = testCase.wvt.variableWithName(['pi_' suffix]);
                testCase.verifyEqual(pressure,testCase.wvt.rho0*testCase.wvt.g*pressureHeight,AbsTol=1e-10);
            end

            total = testCase.wvt.p;
            geostrophic = testCase.wvt.p_g;
            mda = testCase.wvt.p_mda;
            wave = testCase.wvt.p_w;
            inertial = testCase.wvt.p_io;
            testCase.verifyEqual(total,geostrophic+mda+wave+inertial,AbsTol=1e-10);
        end

        function testPeriodicSurfaceInterpolation(testCase)
            field = testCase.wvt.ssh_balanced;
            dx = testCase.wvt.x(2)-testCase.wvt.x(1);
            dy = testCase.wvt.y(2)-testCase.wvt.y(1);
            x = [-0.25*dx 0.4*dx testCase.wvt.Lx-0.2*dx testCase.wvt.Lx+0.3*dx];
            y = [testCase.wvt.Ly-0.3*dy -0.2*dy 0.6*dy testCase.wvt.Ly+0.25*dy];

            xPeriodic = [testCase.wvt.x(end)-testCase.wvt.Lx;testCase.wvt.x;testCase.wvt.Lx];
            yPeriodic = [testCase.wvt.y(end)-testCase.wvt.Ly;testCase.wvt.y;testCase.wvt.Ly];
            fieldPeriodic = field([end 1:end 1],[end 1:end 1]);
            expected = interpn(xPeriodic,yPeriodic,fieldPeriodic,mod(x,testCase.wvt.Lx),mod(y,testCase.wvt.Ly),'linear');
            actual = testCase.wvt.variableAtPositionWithName(x,y,[],'ssh_balanced');
            testCase.verifyEqual(actual,expected,AbsTol=1e-12);
        end

        function testThreeDimensionalInterpolationIsUnchanged(testCase)
            ix = [1 3 testCase.wvt.Nx];
            iy = [testCase.wvt.Ny 4 1];
            iz = [1 3 testCase.wvt.Nz];
            x = reshape(testCase.wvt.x(ix),1,[]) + [-testCase.wvt.Lx 0 testCase.wvt.Lx];
            y = reshape(testCase.wvt.y(iy),1,[]) + [testCase.wvt.Ly 0 -testCase.wvt.Ly];
            z = reshape(testCase.wvt.z(iz),1,[]);

            field = testCase.wvt.u_balanced;
            expected = field(sub2ind(size(field),ix,iy,iz));
            actual = testCase.wvt.variableAtPositionWithName(x,y,z,'u_balanced');
            testCase.verifyEqual(actual,expected,AbsTol=1e-12);
        end
    end

end
