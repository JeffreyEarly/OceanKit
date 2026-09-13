classdef TestWVCompiledConstantPrimitives < matlab.unittest.TestCase
    % Parity tests for constant-stratification compiled primitives.
    methods (Test,TestTags="optional")
        function rawVerticalOperationsMatchMatlab(testCase)
            for isHydrostatic = [true false]
                wvt = makeTransform(isHydrostatic); cleanup = onCleanup(@()delete(wvt));
                backend = WVCompiledTransformBackend.create(wvt); backendCleanup = onCleanup(@()delete(backend));
                modal = complexDeterministic([wvt.Nj wvt.Nkl]);
                grid = complexDeterministic([wvt.Nz wvt.Nkl]);
                Fw = wvt.F_g./wvt.F_wg;
                Gw = wvt.G_g./wvt.G_wg;
                operations = { ...
                    'reconstructF',modal,wvt.iDCT*(wvt.F_g.*modal); ...
                    'projectF',grid,(wvt.DCT*grid)./wvt.F_g; ...
                    'reconstructG',modal,wvt.iDST*(wvt.G_g.*modal); ...
                    'projectG',grid,(wvt.DST*grid)./wvt.G_g; ...
                    'reconstructFw',modal,wvt.iDCT*(Fw.*modal); ...
                    'projectFw',grid,(wvt.DCT*grid)./Fw; ...
                    'reconstructGw',modal,wvt.iDST*(Gw.*modal); ...
                    'projectGw',grid,(wvt.DST*grid)./Gw};
                for index = 1:size(operations,1)
                    options = struct("operator",char(operations{index,1}));
                    result = backend.operation(wvt,"vertical",{operations{index,2}},options);
                    verifyOperationResult(testCase,result,operations{index,3});
                end
                options = struct("operator",'balancedGToWaveG');
                verifyOperationResult(testCase,backend.operation(wvt,"vertical",{modal},options),wvt.transformWithG_wg(modal));

                zeroColumn = complex(zeros(wvt.Nz,1));
                options = struct("operator",'projectFw',"column",1);
                result = backend.operation(wvt,"verticalColumn",{zeroColumn},options);
                verifyOperationResult(testCase,result,wvt.transformFromSpatialDomainWithFio(zeroColumn));
                testCase.verifyEqual(result.values{1},complex(zeros(wvt.Nj,1)));
                clear backendCleanup cleanup
            end
        end

        function horizontalOperationsMatchMatlab(testCase)
            for isHydrostatic = [true false]
                wvt = makeTransform(isHydrostatic); cleanup = onCleanup(@()delete(wvt));
                backend = WVCompiledTransformBackend.create(wvt); backendCleanup = onCleanup(@()delete(backend));
                volume = horizontalFixture(wvt);
                spectrum = wvt.transformFromSpatialDomainWithFourier(volume);
                verifyOperationResult(testCase,backend.operation(wvt,"horizontalForward",{volume}),spectrum);
                expected = wvt.transformToSpatialDomainWithFourier(spectrum);
                verifyOperationResult(testCase,backend.operation(wvt,"horizontalInverse",{complex(spectrum)}),expected);
                for direction = ["x" "y"]
                    for order = [1 2 4]
                        if direction == "x"
                            expected = wvt.diffX(volume,n=order);
                        else
                            expected = wvt.diffY(volume,n=order);
                        end
                        options = struct("direction",char(direction),"order",order);
                        verifyOperationResult(testCase,backend.operation(wvt,"differentiateHorizontal",{volume},options),expected);
                    end
                end
                clear backendCleanup cleanup
            end
        end

        function publicProjectionMatchesMatlabPhaseAndZeroMode(testCase)
            for isHydrostatic = [true false]
                wvt = makeTransform(isHydrostatic); cleanup = onCleanup(@()delete(wvt));
                backend = WVCompiledTransformBackend.create(wvt); backendCleanup = onCleanup(@()delete(backend));
                [u,v,w,eta] = projectionFixture(wvt);
                wvt.t0 = 5; wvt.t = 13;
                if isHydrostatic
                    expected = cell(1,3); [expected{:}] = wvt.transformUVEtaToWaveVortex(u,v,eta);
                    inputs = {u,v,eta};
                else
                    expected = cell(1,3); [expected{:}] = wvt.transformUVWEtaToWaveVortex(u,v,w,eta);
                    inputs = {u,v,w,eta};
                end
                result = backend.operation(wvt,"toWaveVortex",inputs,struct("t",13,"t0",5));
                testCase.verifyEqual(numel(result.values),3);
                for index = 1:3
                    verifyNear(testCase,result.values{index},expected{index});
                    verifyNear(testCase,result.values{index}(:,1),expected{index}(:,1));
                end
                testCase.verifyEqual(result.values{2}(:,1),conj(result.values{1}(:,1)),AbsTol=2*eps);
                testCase.verifyGreaterThan(max(abs(result.values{1}(:,1))),0);
                clear backendCleanup cleanup
            end
        end

        function fieldsAndFluxMatchMatlabSeparately(testCase)
            for isHydrostatic = [true false]
                wvt = makeTransform(isHydrostatic); cleanup = onCleanup(@()delete(wvt));
                initializeInteractingState(wvt);
                backend = WVCompiledTransformBackend.create(wvt); backendCleanup = onCleanup(@()delete(backend));
                names = {'u','v','w','eta'};
                backend.prepare(names);
                reference = cell(size(names)); [reference{:}] = wvt.variableWithName(names{:});
                fields = backend.evaluate(wvt,names);
                testCase.verifyEmpty(fields.flux);
                for index = 1:numel(names)
                    verifyNear(testCase,fields.values{index},reference{index});
                end

                fluxResult = backend.evaluate(wvt,names,shouldEvaluateFlux=true);
                expectedFlux = cell(1,3); [expectedFlux{:}] = wvt.nonlinearFlux();
                testCase.verifyEqual(numel(fluxResult.flux),3);
                fluxMagnitude = 0;
                for index = 1:3
                    verifyNear(testCase,fluxResult.flux{index},expectedFlux{index});
                    fluxMagnitude = max(fluxMagnitude,max(abs(expectedFlux{index}(:))));
                end
                testCase.verifyGreaterThan(fluxMagnitude,0);
                for index = 1:numel(names)
                    verifyNear(testCase,fluxResult.values{index},reference{index});
                end
                clear backendCleanup cleanup
            end
        end
    end
end

function wvt = makeTransform(isHydrostatic)
wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 7],N0=5.2e-3,rho0=1027,latitude=33,isHydrostatic=isHydrostatic,shouldAntialias=true);
end

function value = complexDeterministic(shape)
rng(505,"twister");
value = complex(rand(shape),rand(shape));
end

function volume = horizontalFixture(wvt)
[x,y,z] = ndgrid(0:wvt.Nx-1,0:wvt.Ny-1,0:wvt.Nz-1);
volume = 2+cos(2*pi*x/wvt.Nx)+0.2*sin(4*pi*y/wvt.Ny)+0.1*cos(pi*x)+0.13*cos(pi*y)+0.07*cos(2*pi*(x/wvt.Nx+y/wvt.Ny))+0.01*z;
end

function [u,v,w,eta] = projectionFixture(wvt)
[x,y,z] = ndgrid(0:wvt.Nx-1,0:wvt.Ny-1,0:wvt.Nz-1);
vertical = cos(pi*z/(wvt.Nz-1));
u = 0.4+0.03*vertical+0.2*cos(2*pi*x/wvt.Nx)+0.07*sin(2*pi*(x/wvt.Nx+y/wvt.Ny));
v = -0.2+0.05*vertical+0.13*sin(2*pi*y/wvt.Ny)+0.04*cos(2*pi*(x/wvt.Nx-y/wvt.Ny));
w = 0.03*sin(pi*z/(wvt.Nz-1)).*cos(2*pi*(x/wvt.Nx+y/wvt.Ny));
eta = 0.6*cos(pi*z/(wvt.Nz-1))+0.09*sin(2*pi*x/wvt.Nx).*sin(pi*z/(wvt.Nz-1));
end

function initializeInteractingState(wvt)
wvt.Ap(:) = 0; wvt.Am(:) = 0; wvt.A0(:) = 0;
first = wvt.indexFromModeNumber(1,0,1);
second = wvt.indexFromModeNumber(0,1,2);
third = wvt.indexFromModeNumber(1,1,1);
wvt.Ap(first) = 2.0e-3+0.7e-3i;
wvt.Am(second) = -1.3e-3+0.4e-3i;
wvt.A0(third) = 1.1e-3-0.5e-3i;
wvt.Ap(second) = 0.6e-3-0.3e-3i;
wvt.A0(first) = -0.8e-3+0.2e-3i;
wvt.t0 = 5; wvt.t = 13;
end

function verifyOperationResult(testCase,result,expected)
testCase.verifyTrue(isstruct(result));
testCase.verifyTrue(iscell(result.values));
testCase.verifyEqual(numel(result.values),1);
verifyNear(testCase,result.values{1},expected);
testCase.verifyTrue(isfield(result,"metrics"));
end

function verifyNear(testCase,actual,expected)
testCase.verifySize(actual,size(expected));
testCase.assertTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))));
error = norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
testCase.verifyLessThanOrEqual(error,1e-12);
end
