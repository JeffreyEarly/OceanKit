classdef TestWVCompiledBarotropicPrimitives < matlab.unittest.TestCase
    methods (Test,TestTags="optional")
        function rawAndStateOperationsMatchMatlab(testCase)
            for mode = [0 1]
                for antialias = [false true]
                    wvt = WVTransformBarotropicQG([4000 3000],[8 6],h=.8,j=mode,shouldAntialias=antialias);
                    cleanup = onCleanup(@()delete(wvt));
                    backend = WVCompiledTransformBackend.create(wvt); bcleanup = onCleanup(@()delete(backend));
                    [x,y] = ndgrid(0:wvt.Nx-1,0:wvt.Ny-1);
                    q = 2+cos(2*pi*x/wvt.Nx)+.2*sin(4*pi*y/wvt.Ny)+.1*cos(pi*x)+.13*sin(2*pi*y/wvt.Ny)+.07*cos(2*pi*(x/wvt.Nx+y/wvt.Ny));
                    expected = wvt.transformFromSpatialDomainWithFourier(q);
                    output = backend.operation(wvt,"horizontalForward",{q}); verifyNear(testCase,output.values{1},expected);
                    output = backend.operation(wvt,"qgPVToA0",{q}); verifyNear(testCase,output.values{1},wvt.transformQGPVToWaveVortex(q));
                    output = backend.operation(wvt,"horizontalInverse",{complex(expected)});
                    verifyNear(testCase,output.values{1},wvt.transformToSpatialDomainWithFourier(expected));
                    for direction = ["x" "y"]
                        for order = [1 2 4]
                            if direction == "x", derivative = wvt.diffX(q,n=order); else, derivative = wvt.diffY(q,n=order); end
                            output = backend.operation(wvt,"differentiateHorizontal",{q},struct("direction",char(direction),"order",order));
                            verifyNear(testCase,output.values{1},derivative);
                        end
                    end
                    wvt.A0 = 1e-5*expected; wvt.t = 13;
                    backend.prepare({'u','v'}); mixed = backend.evaluate(wvt,{'u','v'},shouldEvaluateFlux=true);
                    verifyNear(testCase,mixed.values{1},wvt.u); verifyNear(testCase,mixed.values{2},wvt.v);
                    testCase.verifyEqual(numel(mixed.flux),1); verifyNear(testCase,mixed.flux{1},wvt.nonlinearFlux());
                    before = backend.metadata(); fields = backend.evaluate(wvt,{'u','v'});
                    verifyNear(testCase,fields.values{1},wvt.u); verifyNear(testCase,fields.values{2},wvt.v);
                    after = backend.metadata();
                    testCase.verifyEqual(after.engineBytes,before.engineBytes);
                    testCase.verifyEqual(after.sourceBytes,0); testCase.verifyEqual(after.phasePreparations,0);
                    testCase.verifyEqual(after.duplicateExecutions,0); testCase.verifyEqual(after.liveEvaluationBytes,0);
                    oldDepth = wvt.h;
                    testCase.verifyError(@()assignDepth(wvt,2*oldDepth),"WaveVortexModel:CompiledTransformConfigurationLocked");
                    testCase.verifyEqual(wvt.h,oldDepth);
                    restored = backend.evaluate(wvt,{'u','v'});
                    verifyNear(testCase,restored.values{1},fields.values{1});
                    clear bcleanup cleanup
                end
            end
        end
    end
end

function assignDepth(wvt,value)
wvt.h = value;
end

function verifyNear(testCase,actual,expected)
testCase.verifySize(actual,size(expected));
testCase.assertTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))));
error = norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
testCase.verifyLessThanOrEqual(error,1e-12);
end
