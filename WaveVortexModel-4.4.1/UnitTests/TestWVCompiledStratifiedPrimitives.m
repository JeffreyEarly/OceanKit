classdef TestWVCompiledStratifiedPrimitives < matlab.unittest.TestCase
    % Parity tests for the private raw-operation ABI.
    methods (Test,TestTags="optional")
        function stratifiedFieldsAndFlux(testCase)
            for family = ["H" "B" "SQG"]
                wvt = makeTransform(family); cleanup = onCleanup(@()delete(wvt));
                rng(504,"twister"); wvt.initWithRandomFlow(uvMax=.01); wvt.t0 = 5; wvt.t = 13;
                backend = WVCompiledTransformBackend.create(wvt); bcleanup = onCleanup(@()delete(backend));
                names = {'u','v','w','eta'}; flux = cell(1,3);
                if family == "SQG", names = {'u','v'}; flux = cell(1,1); end
                backend.prepare(names);
                actual = backend.evaluate(wvt,names,shouldEvaluateFlux=true);
                reference = cell(size(names)); [reference{:}] = wvt.variableWithName(names{:});
                [flux{:}] = wvt.nonlinearFlux();
                for k = 1:numel(names), testCase.verifyLessThan(relativeError(actual.values{k},reference{k}),1e-12); end
                for k = 1:numel(flux), testCase.verifyLessThan(relativeError(actual.flux{k},flux{k}),1e-12); end
                before = backend.metadata(); fields = backend.evaluate(wvt,names);
                for k = 1:numel(names), testCase.verifyLessThan(relativeError(fields.values{k},reference{k}),1e-12); end
                after = backend.metadata();
                testCase.verifyEqual(after.engineBytes,before.engineBytes);
                testCase.verifyEqual(after.duplicateExecutions,0); testCase.verifyEqual(after.liveEvaluationBytes,0);
                clear bcleanup cleanup
            end
        end

        function hydrostaticVerticalAndLifecycle(testCase)
            wvt = makeTransform("H"); cleanup = onCleanup(@()delete(wvt));
            backend = WVCompiledTransformBackend.create(wvt); bcleanup = onCleanup(@()delete(backend));
            modal = complexDeterministic([wvt.Nj wvt.Nkl]); grid = complexDeterministic([wvt.Nz wvt.Nkl]);
            for operator = ["reconstructF" "projectF" "reconstructG" "projectG"]
                if startsWith(operator,"reconstruct")
                    input = modal;
                    if operator == "reconstructF", expected = wvt.PF0inv*(wvt.P0.*modal); else, expected = wvt.QG0inv*(wvt.Q0.*modal); end
                else
                    input = grid;
                    if operator == "projectF", expected = (wvt.PF0*grid)./wvt.P0; else, expected = (wvt.QG0*grid)./wvt.Q0; end
                end
                result = backend.operation(wvt,"vertical",{input},struct("operator",char(operator)));
                verifyOperationResult(testCase,result,expected);
            end
            before = wv_compiled_backend_mex('moduleMetrics');
            other = makeTransform("H"); otherCleanup = onCleanup(@()delete(other));
            testCase.verifyError(@()backend.operation(other,"vertical",{modal},struct("operator",'reconstructF')),"WaveVortexModel:CompiledTransformMismatch");
            testCase.verifyError(@()backend.operation(wvt,"vertical",{ones(wvt.Nj+1,wvt.Nkl)},struct("operator",'reconstructF')),"WaveVortexModel:CompiledTransformInput");
            testCase.verifyError(@()backend.operation(wvt,"vertical",{modal},struct("operator",'bogus')),"WaveVortexModel:CompiledTransformInput");
            after = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(after.activePlans,before.activePlans);
            clear bcleanup cleanup
        end

        function boussinesqWavePrimitives(testCase)
            wvt = makeTransform("B"); cleanup = onCleanup(@()delete(wvt));
            backend = WVCompiledTransformBackend.create(wvt); bcleanup = onCleanup(@()delete(backend));
            modal = complexDeterministic([wvt.Nj wvt.Nkl]);
            for operator = ["reconstructFw" "reconstructGw"]
                input = modal;
                if operator == "reconstructFw", expected = wvt.transformToSpatialDomainWithFw(modal); else, expected = wvt.transformToSpatialDomainWithGw(modal); end
                verifyOperationResult(testCase,backend.operation(wvt,"vertical",{input},struct("operator",char(operator))),expected);
            end
            result = backend.operation(wvt,"vertical",{modal},struct("operator",'balancedGToWaveG'));
            verifyOperationResult(testCase,result,wvt.transformWithG_wg(modal));
            clear bcleanup cleanup
        end

        function verticalColumnsAndHorizontalOrders(testCase)
            wvt = makeTransform("H"); cleanup = onCleanup(@()delete(wvt));
            backend = WVCompiledTransformBackend.create(wvt); bcleanup = onCleanup(@()delete(backend));
            modal = complexDeterministic([wvt.Nj wvt.Nkl]);
            modal(:,1) = 0;
            for column = [1 wvt.Nkl]
                expected = wvt.PF0inv*(wvt.P0(:).*modal(:,column));
                result = backend.operation(wvt,"verticalColumn",{complex(modal(:,column))},struct("operator",'reconstructF',"column",column));
                verifyOperationResult(testCase,result,expected);
                if column == 1, testCase.verifyEqual(result.values{1},complex(zeros(wvt.Nz,1))); end
            end
            volume = real(complexDeterministic(wvt.spatialMatrixSize));
            for direction = ["x" "y"]
                for order = [1 2 4]
                    if direction == "x", expected = wvt.diffX(volume,n=order); else, expected = wvt.diffY(volume,n=order); end
                    result = backend.operation(wvt,"differentiateHorizontal",{volume},struct("direction",char(direction),"order",order));
                    verifyOperationResult(testCase,result,expected);
                end
            end
            spectrum = wvt.transformFromSpatialDomainWithFourier(volume);
            verifyOperationResult(testCase,backend.operation(wvt,"horizontalForward",{volume},struct()),spectrum);
            verifyOperationResult(testCase,backend.operation(wvt,"horizontalInverse",{spectrum},struct()),wvt.transformToSpatialDomainWithFourier(spectrum));
            clear bcleanup cleanup
        end

        function stateAndQGPVOperations(testCase)
            for family = ["H" "B" "SQG"]
                wvt = makeTransform(family); cleanup = onCleanup(@()delete(wvt));
                backend = WVCompiledTransformBackend.create(wvt); bcleanup = onCleanup(@()delete(backend));
                u = real(complexDeterministic(wvt.spatialMatrixSize)); v = circshift(u,[1 0 0]); eta = 0.1*u;
                wvt.t0 = 5; wvt.t = 13;
                if family == "B"
                    w = 0.1*circshift(u,[0 0 1]); expected = cell(1,3); [expected{:}] = wvt.transformUVWEtaToWaveVortex(u,v,w,eta); inputs = {u,v,w,eta};
                elseif family == "H"
                    expected = cell(1,3); [expected{:}] = wvt.transformUVEtaToWaveVortex(u,v,eta); inputs = {u,v,eta};
                else
                    [~,~,a0] = wvt.transformUVEtaToWaveVortex(u,v,eta); expected = {a0}; inputs = {u,v,eta};
                end
                result = backend.operation(wvt,"toWaveVortex",inputs,struct("t",13,"t0",5));
                testCase.verifyEqual(numel(result.values),numel(expected));
                for i = 1:numel(expected), testCase.verifyLessThanOrEqual(relativeError(result.values{i},expected{i}),1e-12); end
                if family == "SQG"
                    qgpv = real(complexDeterministic(wvt.spatialMatrixSize)); expectedA0 = wvt.transformQGPVToWaveVortex(qgpv);
                    verifyOperationResult(testCase,backend.operation(wvt,"qgPVToA0",{qgpv},struct()),expectedA0);
                end
                testCase.verifyError(@()backend.operation(wvt,"toWaveVortex",{u},struct("t",13,"t0",5)),"WaveVortexModel:CompiledTransformInput");
                clear bcleanup cleanup
            end
        end
    end
end

function wvt = makeTransform(family)
profile = @(z) 2e-5*exp(2*z/1300);
switch family
    case "H", wvt = WVTransformHydrostatic([4000 3000 1000],[8 6 6],Nj=3,N2Function=profile,shouldAntialias=false);
    case "B", wvt = WVTransformBoussinesq([4000 3000 1000],[8 6 6],Nj=3,N2Function=profile,shouldAntialias=false);
    otherwise, wvt = WVTransformStratifiedQG([4000 3000 1000],[8 6 6],Nj=3,N2Function=profile,shouldAntialias=false);
end
end

function value = complexDeterministic(shape)
rng(4001,"twister"); value = rand(shape)+1i*rand(shape);
end

function verifyOperationResult(testCase,result,expected)
testCase.verifyTrue(isstruct(result)); testCase.verifyTrue(iscell(result.values));
testCase.verifyEqual(numel(result.values),1);
testCase.verifyLessThanOrEqual(relativeError(result.values{1},expected),1e-12);
testCase.verifyTrue(isfield(result,"metrics"));
end

function value = relativeError(actual,expected)
if ~all(isfinite(actual(:))) || ~all(isfinite(expected(:)))
    value = Inf;
else
    value = max(abs(actual(:)-expected(:)))/max(max(abs(expected(:))),realmin);
end
end
