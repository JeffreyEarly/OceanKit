classdef TestWVFourierStorageLayout < matlab.unittest.TestCase
    methods (Test,TestTags="full")
        function fullHalfLayoutsAreEquivalent(testCase)
            sizes = [8 6 3;9 7 4];
            for iSize = 1:size(sizes,1)
                Nx = sizes(iSize,1);
                Ny = sizes(iSize,2);
                nBatch = sizes(iSize,3);
                for conjugateDimension = [1 2]
                    for shouldAntialias = [false true]
                        for shouldExcludeNyquist = [false true]
                            geometry = WVGeometryDoublyPeriodic([4000 3000],[Nx Ny],Nz=nBatch,conjugateDimension=conjugateDimension,shouldAntialias=shouldAntialias,shouldExcludeNyquist=shouldExcludeNyquist,shouldExcludeConjugates=true);
                            fullLayout = WVFourierStorageLayout(geometry,"full-complex");
                            halfXLayout = WVFourierStorageLayout(geometry,"hermitian-half",compressedDimension=1);
                            halfYLayout = WVFourierStorageLayout(geometry,"hermitian-half",compressedDimension=2);

                            rng(7000+100*iSize+10*conjugateDimension+2*shouldAntialias+shouldExcludeNyquist,"twister");
                            realInput = randn(Nx,Ny,nBatch);
                            fullSpectrum = fft(fft(realInput,Nx,1),Ny,2)/(Nx*Ny);
                            canonical = fullLayout.transformFromFourierStorageToWVGrid(fullLayout.reshapeFourierStorageToRows(fullSpectrum));
                            halfXSpectrum = fullSpectrum(1:(floor(Nx/2)+1),:,:);
                            halfYSpectrum = fullSpectrum(:,1:(floor(Ny/2)+1),:);

                            testCase.verifyEqual(halfXLayout.transformFromFourierStorageToWVGrid(halfXLayout.reshapeFourierStorageToRows(halfXSpectrum)),canonical,AbsTol=1e-12);
                            testCase.verifyEqual(halfYLayout.transformFromFourierStorageToWVGrid(halfYLayout.reshapeFourierStorageToRows(halfYSpectrum)),canonical,AbsTol=1e-12);

                            fullRows = fullLayout.transformFromWVGridToFourierStorage(fullLayout.allocateFourierStorage(nBatch),canonical);
                            halfXRows = halfXLayout.transformFromWVGridToFourierStorage(halfXLayout.allocateFourierStorage(nBatch),canonical);
                            halfYRows = halfYLayout.transformFromWVGridToFourierStorage(halfYLayout.allocateFourierStorage(nBatch),canonical);
                            testCase.verifyEqual(fullLayout.transformFromFourierStorageToWVGrid(fullRows),canonical,AbsTol=1e-12);
                            testCase.verifyEqual(halfXLayout.transformFromFourierStorageToWVGrid(halfXRows),canonical,AbsTol=1e-12);
                            testCase.verifyEqual(halfYLayout.transformFromFourierStorageToWVGrid(halfYRows),canonical,AbsTol=1e-12);

                            fullReal = inverseFromFullLayout(fullLayout,fullRows,conjugateDimension,Nx,Ny);
                            halfXReal = inverseFromHalfLayout(halfXLayout,halfXRows,Nx,Ny);
                            halfYReal = inverseFromHalfLayout(halfYLayout,halfYRows,Nx,Ny);
                            testCase.verifyEqual(halfXReal,fullReal,AbsTol=1e-12);
                            testCase.verifyEqual(halfYReal,fullReal,AbsTol=1e-12);
                        end
                    end
                end
            end
        end

        function layoutMetadataAndStorageAreExact(testCase)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[10 7],Nz=5,shouldAntialias=false,shouldExcludeNyquist=false);
            layouts = [WVFourierStorageLayout(geometry,"full-complex") WVFourierStorageLayout(geometry,"hermitian-half",compressedDimension=1) WVFourierStorageLayout(geometry,"hermitian-half",compressedDimension=2)];
            expectedShapes = [10 7;6 7;10 4];
            for iLayout = 1:numel(layouts)
                layout = layouts(iLayout);
                testCase.verifyEqual(layout.fourierStorageSize,expectedShapes(iLayout,:));
                testCase.verifyEqual(layout.nFourierStorageRows,prod(expectedShapes(iLayout,:)));
                testCase.verifyEqual(layout.mappingMethod,"two-dimensional-rows");
                ledger = layout.mappingMemoryUsage();
                testCase.verifyEqual(layout.mappingMemoryBytes,sum([ledger.bytes]));
                testCase.verifyTrue(all(string({ledger.class}) == "uint64"));
                for entry = ledger
                    testCase.verifyEqual(entry.bytes,8*prod(entry.shape));
                end
                rows = layout.allocateFourierStorage(5);
                testCase.verifySize(rows,[layout.nFourierStorageRows 5]);
                storage = layout.reshapeFourierRowsToStorage(rows);
                testCase.verifyEqual(size(storage),[layout.fourierStorageSize 5]);
                testCase.verifyEqual(layout.reshapeFourierStorageToRows(storage),rows);
                testCase.verifyTrue(all([layout.fourierRowsForDirectWVIndices;layout.fourierRowsForConjugatedWVIndices] >= 1));
                testCase.verifyTrue(all([layout.fourierRowsForDirectWVIndices;layout.fourierRowsForConjugatedWVIndices] <= layout.nFourierStorageRows));
            end
        end

        function selfConjugateEntriesAreRealOnInsertion(testCase)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[8 6],Nz=2,shouldAntialias=false,shouldExcludeNyquist=false);
            for compressedDimension = [1 2]
                layout = WVFourierStorageLayout(geometry,"hermitian-half",compressedDimension=compressedDimension);
                canonical = complex(randn(2,geometry.Nkl),randn(2,geometry.Nkl));
                rows = layout.transformFromWVGridToFourierStorage(layout.allocateFourierStorage(2),canonical);
                testCase.verifyEqual(imag(rows(layout.selfConjugateFourierRows,:)),zeros(numel(layout.selfConjugateFourierRows),2));
            end
        end

        function invalidShapesAndConfigurationsAreRejected(testCase)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[8 6],Nz=2);
            testCase.verifyError(@()WVFourierStorageLayout(geometry,"full-complex",compressedDimension=1),"WaveVortexModel:InvalidFullFourierStorageLayout");
            testCase.verifyError(@()WVFourierStorageLayout(geometry,"hermitian-half"),"WaveVortexModel:InvalidHalfFourierStorageLayout");
            layout = WVFourierStorageLayout(geometry,"hermitian-half",compressedDimension=1);
            testCase.verifyError(@()layout.transformFromFourierStorageToWVGrid(complex(zeros(layout.nFourierStorageRows+1,2))),"WaveVortexModel:InvalidFourierStorageShape");
            testCase.verifyError(@()layout.transformFromWVGridToFourierStorage(layout.allocateFourierStorage(2),complex(zeros(3,geometry.Nkl))),"WaveVortexModel:InvalidWVGridShape");
            testCase.verifyError(@()layout.reshapeFourierStorageToRows(complex(zeros(8,6,2))),"WaveVortexModel:InvalidFourierStorageShape");
        end

        function classIsVisibleSealedReadOnlyAndDocumented(testCase)
            classMetadata = meta.class.fromName("WVFourierStorageLayout");
            testCase.verifyNotEmpty(classMetadata);
            testCase.verifyFalse(classMetadata.Hidden);
            testCase.verifyTrue(classMetadata.Sealed);
            testCase.verifyNotEmpty(classMetadata.DetailedDescription);

            expectedProperties = ["Nkl" "compressedDimension" "conjugatedWVIndices" ...
                "directWVIndices" "fourierRowsForConjugatedWVIndices" ...
                "fourierRowsForDirectWVIndices" "fourierStorageSize" ...
                "fourierStorageType" "hermitianCompletionRows" ...
                "hermitianSourceRows" "hermitianSourceWVIndices" ...
                "horizontalGridSize" "mappingMemoryBytes" "mappingMethod" ...
                "nFourierStorageRows" "selfConjugateFourierRows"];
            properties = classMetadata.PropertyList;
            if iscell(properties)
                properties = [properties{:}];
            end
            publicProperties = properties(strcmp({properties.GetAccess},"public"));
            testCase.verifyEqual(sort(string({publicProperties.Name})),sort(expectedProperties));
            testCase.verifyTrue(all(strcmp({publicProperties.SetAccess},"private")));

            expectedMethods = ["WVFourierStorageLayout" "allocateFourierStorage" ...
                "expandedLegacyMappings" "mappingMemoryUsage" ...
                "reshapeFourierRowsToStorage" "reshapeFourierStorageToRows" ...
                "transformFromFourierStorageToWVGrid" ...
                "transformFromWVGridToFourierStorage"];
            methods = classMetadata.MethodList;
            if iscell(methods)
                methods = [methods{:}];
            end
            publicMethods = methods(strcmp({methods.Access},"public") & ~strcmp({methods.Name},"empty"));
            testCase.verifyEqual(sort(string({publicMethods.Name})),sort(expectedMethods));
            testCase.verifyTrue(all(strlength(string({publicMethods.DetailedDescription})) > 0));
        end

        function callerOwnedStorageMustBeReassigned(testCase)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[8 6],Nz=2,shouldAntialias=false,shouldExcludeNyquist=false);
            layout = WVFourierStorageLayout(geometry,"full-complex");
            originalRows = layout.allocateFourierStorage(2);
            rows = originalRows;
            wvArray = complex(randn(2,geometry.Nkl),randn(2,geometry.Nkl));
            layout.transformFromWVGridToFourierStorage(rows,wvArray);
            testCase.verifyEqual(rows,originalRows);
            rows = layout.transformFromWVGridToFourierStorage(rows,wvArray);
            testCase.verifyNotEqual(rows,originalRows);
            testCase.verifyEqual(layout.transformFromFourierStorageToWVGrid(rows),wvArray,AbsTol=1e-12);
        end

        function builtinUsesRowsAndPreservesLegacyMappings(testCase)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[16 12],Nz=5,shouldAntialias=false);
            diagnostics = geometry.fourierStorageLayoutDiagnostics();
            testCase.verifyEqual(diagnostics.fourierStorageType,"full-complex");
            testCase.verifyEqual(diagnostics.mappingMethod,"two-dimensional-rows");
            testCase.verifyFalse(diagnostics.legacyMappingsAreMaterialized);
            testCase.verifyEqual(diagnostics.legacyMappingBytes,0);

            rng(7100,"twister");
            realInput = randn(16,12,5);
            canonical = geometry.transformFromSpatialDomainWithFourier(realInput);
            realOutput = geometry.transformToSpatialDomainWithFourier(canonical);
            testCase.verifySize(canonical,[5 geometry.Nkl]);
            testCase.verifySize(realOutput,[16 12 5]);
            testCase.verifySize(geometry.fastTransform.complexBuffer,[16 12 5]);
            diagnostics = geometry.fourierStorageLayoutDiagnostics();
            testCase.verifyFalse(diagnostics.legacyMappingsAreMaterialized);

            [expectedPrimary,expectedConjugate,expectedWVConjugate] = geometry.indicesFromWVGridToDFTGrid(5,isHalfComplex=true);
            testCase.verifyEqual(geometry.dftPrimaryIndex,uint64(expectedPrimary));
            testCase.verifyEqual(geometry.dftConjugateIndex,uint64(expectedConjugate));
            testCase.verifyEqual(geometry.wvConjugateIndex,uint64(expectedWVConjugate));
            diagnostics = geometry.fourierStorageLayoutDiagnostics();
            testCase.verifyTrue(diagnostics.legacyMappingsAreMaterialized);
            testCase.verifyEqual(diagnostics.legacyMappingBytes,8*(numel(expectedPrimary)+numel(expectedConjugate)+numel(expectedWVConjugate)));
            testCase.verifyEqual(geometry.dftPrimaryIndex,uint64(expectedPrimary));
        end

        function builtinDefaultPathMatchesFrozenExpressions(testCase)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[16 12],Nz=5,shouldAntialias=false);
            rng(7200,"twister");
            realInput = randn(16,12,5);
            fullSpectrum = fft(fft(realInput,16,1),12,2)/(16*12);
            [primary,conjugate,wvConjugate] = geometry.indicesFromWVGridToDFTGrid(5,isHalfComplex=true);
            expectedCanonical = reshape(fullSpectrum(primary),[5 geometry.Nkl]);
            actualCanonical = geometry.transformFromSpatialDomainWithFourier(realInput);
            testCase.verifyEqual(actualCanonical,expectedCanonical,AbsTol=1e-12);

            expectedBuffer = complex(zeros(16,12,5));
            expectedBuffer(primary) = expectedCanonical;
            expectedBuffer(conjugate) = conj(expectedCanonical(wvConjugate));
            expectedReal = ifft(ifft(expectedBuffer,16,1),12,2,"symmetric")*(16*12);
            actualReal = geometry.transformToSpatialDomainWithFourier(actualCanonical);
            testCase.verifyEqual(actualReal,expectedReal,AbsTol=1e-12);
            testCase.verifyEqual(geometry.fastTransform.complexBuffer,expectedBuffer,AbsTol=1e-12);
        end
    end
end

function realOutput = inverseFromFullLayout(layout,rows,conjugateDimension,Nx,Ny)
spectrum = layout.reshapeFourierRowsToStorage(rows);
if conjugateDimension == 1
    realOutput = ifft(ifft(spectrum,Ny,2),Nx,1,"symmetric")*(Nx*Ny);
else
    realOutput = ifft(ifft(spectrum,Nx,1),Ny,2,"symmetric")*(Nx*Ny);
end
end

function realOutput = inverseFromHalfLayout(layout,rows,Nx,Ny)
halfSpectrum = layout.reshapeFourierRowsToStorage(rows);
fullSpectrum = expandHalfSpectrum(halfSpectrum,[Nx Ny],layout.compressedDimension);
realOutput = ifft(ifft(fullSpectrum,Nx,1),Ny,2,"symmetric")*(Nx*Ny);
end

function fullSpectrum = expandHalfSpectrum(halfSpectrum,physicalShape,compressedDimension)
Nx = physicalShape(1);
Ny = physicalShape(2);
nBatch = size(halfSpectrum,3);
fullSpectrum = complex(zeros(Nx,Ny,nBatch));
if compressedDimension == 1
    fullSpectrum(1:size(halfSpectrum,1),:,:) = halfSpectrum;
else
    fullSpectrum(:,1:size(halfSpectrum,2),:) = halfSpectrum;
end
for iK = 1:Nx
    for iL = 1:Ny
        isStored = iK <= size(halfSpectrum,1) && iL <= size(halfSpectrum,2);
        if ~isStored
            conjugateK = mod(Nx-(iK-1),Nx)+1;
            conjugateL = mod(Ny-(iL-1),Ny)+1;
            fullSpectrum(iK,iL,:) = conj(fullSpectrum(conjugateK,conjugateL,:));
        end
    end
end
end
