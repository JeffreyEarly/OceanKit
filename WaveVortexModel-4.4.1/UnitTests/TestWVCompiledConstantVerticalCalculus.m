classdef TestWVCompiledConstantVerticalCalculus < matlab.unittest.TestCase
    % Focused public parity for tiled constant-stratification vertical calculus.
    methods (Test,TestTags="optional")
        function volumesMatchMatlabAcrossTileBoundary(testCase)
            gridSizes = {[8 6 9],[18 16 17]};
            for isHydrostatic = [true false]
                for shouldAntialias = [true false]
                    for iGrid = 1:numel(gridSizes)
                        Nxyz = gridSizes{iGrid};
                        [matlabWVT,compiledWVT] = makeTransforms(Nxyz,isHydrostatic,shouldAntialias);
                        cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                        diagnostic = configurationName(Nxyz,isHydrostatic,shouldAntialias);
                        horizontalRows = prod(Nxyz(1:2));
                        if iGrid == 1
                            testCase.verifyLessThan(horizontalRows,256,diagnostic+" tile-boundary precondition");
                        else
                            testCase.verifyGreaterThan(horizontalRows,256,diagnostic+" tile-boundary precondition");
                        end
                        [F,G,expectedDzF,expectedDzG] = analyticVolumes(matlabWVT);

                        testCase.verifyGreaterThan(norm(expectedDzF(:)),0,diagnostic+" analytic diffZF nonzero");
                        testCase.verifyGreaterThan(norm(expectedDzG(:)),0,diagnostic+" analytic diffZG nonzero");
                        verifyNear(testCase,matlabWVT.diffZF(F),expectedDzF,2e-11,diagnostic+" analytic diffZF");
                        verifyNear(testCase,matlabWVT.diffZG(G),expectedDzG,2e-11,diagnostic+" analytic diffZG");

                        broadband = broadbandVolume(matlabWVT);
                        realF = F+broadband;
                        realG = G-0.3*broadband;
                        valuesF = {realF,complex(realF,0.25*G-0.4*broadband)};
                        valuesG = {realG,complex(realG,-0.25*F+0.2*broadband)};
                        kinds = ["real" "complex"];
                        for iKind = 1:2
                            for order = 1:4
                                verifyDerivative(testCase,matlabWVT,compiledWVT,"diffZF",valuesF{iKind},order,diagnostic+" "+kinds(iKind));
                                verifyDerivative(testCase,matlabWVT,compiledWVT,"diffZG",valuesG{iKind},order,diagnostic+" "+kinds(iKind));
                            end
                            verifyIntegral(testCase,matlabWVT,compiledWVT,"intZF",valuesF{iKind},diagnostic+" "+kinds(iKind)+" volume");
                            verifyIntegral(testCase,matlabWVT,compiledWVT,"intZG",valuesG{iKind},diagnostic+" "+kinds(iKind)+" volume");
                        end
                        clear cleanup
                    end
                end
            end
        end

        function evenGridDerivativesMatchIndependentAnalyticModes(testCase)
            Nxyz = [18 16 16];
            for isHydrostatic = [true false]
                for shouldAntialias = [true false]
                    [matlabWVT,compiledWVT] = makeTransforms(Nxyz,isHydrostatic,shouldAntialias);
                    cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                    diagnostic = configurationName(Nxyz,isHydrostatic,shouldAntialias);
                    [F,G,expectedF,expectedG] = evenGridAnalyticVolumes(matlabWVT);
                    valuesF = {F,complex(F,0.25*F)};
                    valuesG = {G,complex(G,-0.25*G)};
                    kinds = ["real" "complex"];

                    % Mid-spectrum modes keep fourth derivatives well scaled on
                    % the even grid and avoid low-mode cancellation amplification.
                    for iKind = 1:2
                        originalF = valuesF{iKind};
                        originalG = valuesG{iKind};
                        for order = 1:4
                            if iKind == 1
                                expectedFOrder = expectedF{order};
                                expectedGOrder = expectedG{order};
                            else
                                expectedFOrder = complex(expectedF{order},0.25*expectedF{order});
                                expectedGOrder = complex(expectedG{order},-0.25*expectedG{order});
                            end
                            matlabF = matlabWVT.diffZF(valuesF{iKind},n=order);
                            compiledF = compiledWVT.diffZF(valuesF{iKind},n=order);
                            matlabG = matlabWVT.diffZG(valuesG{iKind},n=order);
                            compiledG = compiledWVT.diffZG(valuesG{iKind},n=order);
                            suffix = diagnostic+" "+kinds(iKind)+" order "+order;
                            testCase.verifyGreaterThan(norm(expectedFOrder(:)),0,suffix+" analytic diffZF nonzero");
                            testCase.verifyGreaterThan(norm(expectedGOrder(:)),0,suffix+" analytic diffZG nonzero");
                            verifyNear(testCase,matlabF,expectedFOrder,1e-12,suffix+" MATLAB diffZF analytic");
                            verifyNear(testCase,compiledF,expectedFOrder,1e-12,suffix+" compiled diffZF analytic");
                            verifyNear(testCase,compiledF,matlabF,1e-12,suffix+" diffZF parity");
                            verifyNear(testCase,matlabG,expectedGOrder,1e-12,suffix+" MATLAB diffZG analytic");
                            verifyNear(testCase,compiledG,expectedGOrder,1e-12,suffix+" compiled diffZG analytic");
                            verifyNear(testCase,compiledG,matlabG,1e-12,suffix+" diffZG parity");
                        end
                        testCase.verifyEqual(valuesF{iKind},originalF,diagnostic+" "+kinds(iKind)+" diffZF input unchanged");
                        testCase.verifyEqual(valuesG{iKind},originalG,diagnostic+" "+kinds(iKind)+" diffZG input unchanged");
                    end
                    clear cleanup
                end
            end
        end

        function matrixIntegralsMatchMatlabAtTileTails(testCase)
            columnCounts = [0 1 255 256 257];
            gridSizes = {[8 6 9],[18 16 17],[18 16 16]};
            for isHydrostatic = [true false]
                for shouldAntialias = [true false]
                    for iGrid = 1:numel(gridSizes)
                        Nxyz = gridSizes{iGrid};
                        [matlabWVT,compiledWVT] = makeTransforms(Nxyz,isHydrostatic,shouldAntialias);
                        cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                        diagnostic = configurationName(Nxyz,isHydrostatic,shouldAntialias);
                        for columnCount = columnCounts
                            realMatrix = deterministicMatrix(matlabWVT.Nz,columnCount);
                            complexMatrix = complex(realMatrix,0.5*flip(realMatrix,1));
                            for method = ["intZF" "intZG"]
                                verifyIntegral(testCase,matlabWVT,compiledWVT,method,realMatrix,diagnostic+" real matrix columns="+columnCount);
                                verifyIntegral(testCase,matlabWVT,compiledWVT,method,complexMatrix,diagnostic+" complex matrix columns="+columnCount);
                            end
                        end
                        clear cleanup
                    end
                end
            end
        end

        function tracerEvolutionMatchesMatlab(testCase)
            for isHydrostatic = [true false]
                [matlabWVT,compiledWVT] = makeTransforms([8 6 9],isHydrostatic,true);
                initializeInteractingState(matlabWVT);
                initializeInteractingState(compiledWVT);
                matlabModel = WVModel(matlabWVT);
                compiledModel = WVModel(compiledWVT);
                cleanup = onCleanup(@()deleteModelsAndTransforms(matlabModel,compiledModel,matlabWVT,compiledWVT));
                diagnostic = "tracer evolution, hydrostatic="+string(isHydrostatic);

                [X,Y,Z] = matlabWVT.xyzGrid;
                initialTracer = 1+0.2*sin(2*pi*X/matlabWVT.Lx).*cos(pi*(Z+matlabWVT.Lz)/matlabWVT.Lz);
                initialTracer = initialTracer+0.15*cos(2*pi*Y/matlabWVT.Ly).*sin(2*pi*(Z+matlabWVT.Lz)/matlabWVT.Lz);
                matlabModel.addTracer(initialTracer,"dye");
                compiledModel.addTracer(initialTracer,"dye");
                matlabModel.setupIntegrator(integratorType="fixed",deltaT=1);
                compiledModel.setupIntegrator(integratorType="fixed",deltaT=1);
                matlabModel.integrateToTime(2,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                compiledModel.integrateToTime(2,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);

                expectedTracer = matlabModel.tracer("dye");
                actualTracer = compiledModel.tracer("dye");
                verifyNear(testCase,actualTracer,expectedTracer,1e-12,diagnostic+" tracer");
                verifyNear(testCase,compiledWVT.Ap,matlabWVT.Ap,1e-12,diagnostic+" Ap");
                verifyNear(testCase,compiledWVT.Am,matlabWVT.Am,1e-12,diagnostic+" Am");
                verifyNear(testCase,compiledWVT.A0,matlabWVT.A0,1e-12,diagnostic+" A0");
                testCase.verifyGreaterThan(norm(expectedTracer(:)-initialTracer(:)),1e-12*norm(initialTracer(:)),diagnostic+" nontrivial tracer change");
                clear cleanup
            end
        end
    end
end

function [matlabWVT,compiledWVT] = makeTransforms(Nxyz,isHydrostatic,shouldAntialias)
matlabWVT = WVTransformConstantStratification([6000 4000 1000],Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=isHydrostatic,shouldAntialias=shouldAntialias,computationalBackend="matlab");
compiledWVT = WVTransformConstantStratification([6000 4000 1000],Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=isHydrostatic,shouldAntialias=shouldAntialias,computationalBackend="compiled");
end

function [F,G,expectedDzF,expectedDzG] = analyticVolumes(wvt)
[X,Y,Z] = wvt.xyzGrid;
horizontal = 1+0.2*cos(2*pi*X/wvt.Lx)+0.1*sin(2*pi*Y/wvt.Ly);
mode = 2;
m = mode*pi/wvt.Lz;
theta = m*(Z+wvt.Lz);
F = horizontal.*cos(theta);
G = horizontal.*sin(theta);
expectedDzF = -m*horizontal.*sin(theta);
expectedDzG = m*horizontal.*cos(theta);
end

function [F,G,expectedF,expectedG] = evenGridAnalyticVolumes(wvt)
[X,Y,Z] = wvt.xyzGrid;
horizontal = 1+0.2*cos(2*pi*X/wvt.Lx)+0.1*sin(2*pi*Y/wvt.Ly);
m6 = 6*pi/wvt.Lz;
m3 = 3*pi/wvt.Lz;
theta6 = m6*(Z+wvt.Lz);
theta3 = m3*(Z+wvt.Lz);
cosine = cos(theta6)+0.2*cos(theta3);
sine = sin(theta6)+0.2*sin(theta3);
F = horizontal.*cosine;
G = horizontal.*sine;
expectedF = {
    -horizontal.*(m6*sin(theta6)+0.2*m3*sin(theta3)), ...
    -horizontal.*(m6^2*cos(theta6)+0.2*m3^2*cos(theta3)), ...
    horizontal.*(m6^3*sin(theta6)+0.2*m3^3*sin(theta3)), ...
    horizontal.*(m6^4*cos(theta6)+0.2*m3^4*cos(theta3))};
expectedG = {
    horizontal.*(m6*cos(theta6)+0.2*m3*cos(theta3)), ...
    -horizontal.*(m6^2*sin(theta6)+0.2*m3^2*sin(theta3)), ...
    -horizontal.*(m6^3*cos(theta6)+0.2*m3^3*cos(theta3)), ...
    horizontal.*(m6^4*sin(theta6)+0.2*m3^4*sin(theta3))};
end

function value = broadbandVolume(wvt)
[x,y,z] = ndgrid(0:wvt.Nx-1,0:wvt.Ny-1,0:wvt.Nz-1);
value = 0.013*sin(2*pi*(3*x/wvt.Nx+2*y/wvt.Ny)).*cos(3*pi*z/(wvt.Nz-1));
value = value+0.007*cos(2*pi*(x/wvt.Nx-3*y/wvt.Ny)).*sin(4*pi*z/(wvt.Nz-1));
end

function value = deterministicMatrix(rows,columns)
if columns == 0
    value = zeros(rows,0);
    return
end
[row,column] = ndgrid(0:rows-1,0:columns-1);
value = 0.25+cos(pi*row/max(rows-1,1))+0.1*sin(2*pi*(column+1)/max(columns,2));
end

function initializeInteractingState(wvt)
wvt.Ap(:) = 0;
wvt.Am(:) = 0;
wvt.A0(:) = 0;
first = wvt.indexFromModeNumber(1,0,1);
second = wvt.indexFromModeNumber(0,1,2);
third = wvt.indexFromModeNumber(1,1,1);
wvt.Ap(first) = 2.0e-3+0.7e-3i;
wvt.Am(second) = -1.3e-3+0.4e-3i;
wvt.A0(third) = 1.1e-3-0.5e-3i;
wvt.Ap(second) = 0.6e-3-0.3e-3i;
wvt.A0(first) = -0.8e-3+0.2e-3i;
wvt.t0 = 0;
wvt.t = 0;
end

function verifyDerivative(testCase,matlabWVT,compiledWVT,method,input,order,diagnostic)
original = input;
expected = feval(method,matlabWVT,input,n=order);
actual = feval(method,compiledWVT,input,n=order);
testCase.verifyEqual(input,original,diagnostic+" "+method+" order "+order+" input unchanged");
verifyNear(testCase,actual,expected,1e-12,diagnostic+" "+method+" order "+order);
end

function verifyIntegral(testCase,matlabWVT,compiledWVT,method,input,diagnostic)
original = input;
expected = feval(method,matlabWVT,input);
actual = feval(method,compiledWVT,input);
testCase.verifyEqual(input,original,diagnostic+" "+method+" input unchanged");
verifyNear(testCase,actual,expected,1e-12,diagnostic+" "+method);
end

function verifyNear(testCase,actual,expected,tolerance,diagnostic)
testCase.verifySize(actual,size(expected),diagnostic);
testCase.assertTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))),diagnostic+" finite output");
relativeError = norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
testCase.verifyLessThanOrEqual(relativeError,tolerance,diagnostic);
end

function value = configurationName(Nxyz,isHydrostatic,shouldAntialias)
value = "grid="+join(string(Nxyz),"x")+", hydrostatic="+string(isHydrostatic)+", antialias="+string(shouldAntialias);
end

function deleteTransforms(varargin)
for transform = varargin
    if isvalid(transform{1}), delete(transform{1}), end
end
end

function deleteModelsAndTransforms(matlabModel,compiledModel,matlabWVT,compiledWVT)
if isvalid(matlabModel), delete(matlabModel), end
if isvalid(compiledModel), delete(compiledModel), end
deleteTransforms(matlabWVT,compiledWVT);
end
