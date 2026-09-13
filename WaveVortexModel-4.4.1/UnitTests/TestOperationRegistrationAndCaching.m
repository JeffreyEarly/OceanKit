classdef TestOperationRegistrationAndCaching < matlab.unittest.TestCase
    methods (Test, TestTags="full")
        function singleAndMultipleOutputsExecuteOnceAndRemainOrdered(testCase)
            wvt = TestOperationRegistrationAndCaching.transform();
            singleAnnotation = TestOperationRegistrationAndCaching.annotation("single");
            single = WVCountingOperation('single',singleAnnotation,@(~) {11});
            pairAnnotations = [TestOperationRegistrationAndCaching.annotation("pair_a") TestOperationRegistrationAndCaching.annotation("pair_b")];
            pair = WVCountingOperation('pair',pairAnnotations,@(~) {21,22});
            wvt.addOperation([single pair]);

            testCase.verifyTrue(wvt.operationWithName('single') == single)
            testCase.verifyTrue(wvt.operationWithName('pair') == pair)
            testCase.verifyTrue(wvt.propertyAnnotationWithName("single").modelOp == single)
            testCase.verifyTrue(wvt.propertyAnnotationWithName("pair_a").modelOp == pair)

            testCase.verifyEqual(wvt.variableWithName('single'),11)
            testCase.verifyEqual(wvt.single,11)
            testCase.verifyEqual(single.callCount,1)

            [first,second] = wvt.performOperationWithName('pair');
            testCase.verifyEqual([first second],[21 22])
            [firstAgain,secondAgain] = wvt.pair;
            testCase.verifyEqual([firstAgain secondAgain],[21 22])
            [requestedSecond,requestedFirst] = wvt.variableWithName('pair_b','pair_a');
            testCase.verifyEqual([requestedSecond requestedFirst],[22 21])
            testCase.verifyEqual(pair.callCount,1)
            testCase.verifyTrue(all(isKey(wvt.variableCache,["single" "pair_a" "pair_b"])))
        end

        function invalidAndConflictingRegistrationsAreAtomic(testCase)
            wvt = TestOperationRegistrationAndCaching.transform();
            baseline = WVCountingOperation('baseline',TestOperationRegistrationAndCaching.annotation("baseline"),@(~) {31});
            wvt.addOperation(baseline);
            testCase.verifyEqual(wvt.baseline,31)
            initialState = TestOperationRegistrationAndCaching.registryState(wvt);

            unknownDimension = WVVariableAnnotation('unknown_dimension',{'not_a_dimension'},'1','unknown dimension output');
            invalidDimension = WVCountingOperation('unknown_dimension',unknownDimension,@(~) {1});
            testCase.verifyError(@()wvt.addOperation(invalidDimension),'')
            TestOperationRegistrationAndCaching.verifyRegistryState(testCase,wvt,initialState)

            protectedAnnotation = WVCountingOperation('t',TestOperationRegistrationAndCaching.annotation("t"),@(~) {1});
            testCase.verifyError(@()wvt.addOperation(protectedAnnotation,shouldOverwriteExisting=true),'')
            TestOperationRegistrationAndCaching.verifyRegistryState(testCase,wvt,initialState)

            outputConflict = WVCountingOperation('output conflict',[TestOperationRegistrationAndCaching.annotation("baseline") TestOperationRegistrationAndCaching.annotation("conflict_extra")],@(~) {2,3});
            operationNameConflict = WVCountingOperation('baseline',[TestOperationRegistrationAndCaching.annotation("different_output") TestOperationRegistrationAndCaching.annotation("other_output")],@(~) {3,4});
            testCase.verifyError(@()wvt.addOperation(outputConflict),'')
            testCase.verifyError(@()wvt.addOperation(operationNameConflict),'')
            TestOperationRegistrationAndCaching.verifyRegistryState(testCase,wvt,initialState)

            duplicateAnnotations = [TestOperationRegistrationAndCaching.annotation("duplicate") TestOperationRegistrationAndCaching.annotation("duplicate")];
            duplicateOutputs = WVCountingOperation('duplicate outputs',duplicateAnnotations,@(~) {1,2});
            testCase.verifyError(@()wvt.addOperation(duplicateOutputs),'')

            valid = WVCountingOperation('would_be_valid',TestOperationRegistrationAndCaching.annotation("would_be_valid"),@(~) {4});
            testCase.verifyError(@()wvt.addOperation([valid invalidDimension]),'')
            TestOperationRegistrationAndCaching.verifyRegistryState(testCase,wvt,initialState)
            testCase.verifyFalse(ismember("would_be_valid",wvt.variableNames))
        end

        function overwriteRemovesAllConflictsAndPreservesUnrelatedCache(testCase)
            wvt = TestOperationRegistrationAndCaching.transform();
            first = WVCountingOperation('first',[TestOperationRegistrationAndCaching.annotation("first_a") TestOperationRegistrationAndCaching.annotation("shared_a")],@(~) {1,2});
            second = WVCountingOperation('second',[TestOperationRegistrationAndCaching.annotation("second_a") TestOperationRegistrationAndCaching.annotation("shared_b")],@(~) {3,4});
            unrelated = WVCountingOperation('unrelated',TestOperationRegistrationAndCaching.annotation("unrelated"),@(~) {5});
            wvt.addOperation([first second unrelated]);
            wvt.performOperationWithName('first');
            wvt.performOperationWithName('second');
            testCase.verifyEqual(wvt.unrelated,5)

            replacementAnnotations = [TestOperationRegistrationAndCaching.annotation("shared_a") TestOperationRegistrationAndCaching.annotation("shared_b")];
            replacement = WVCountingOperation('replacement',replacementAnnotations,@(~) {12,14});
            wvt.addOperation(replacement,shouldOverwriteExisting=true,shouldSuppressWarning=true);

            testCase.verifyError(@()wvt.operationWithName('first'),'')
            testCase.verifyError(@()wvt.operationWithName('second'),'')
            testCase.verifyFalse(any(ismember(["first_a" "second_a"],wvt.variableNames)))
            testCase.verifyFalse(any(isKey(wvt.variableCache,["first_a" "shared_a" "second_a" "shared_b"])))
            testCase.verifyTrue(isKey(wvt.variableCache,"unrelated"))
            [sharedA,sharedB] = wvt.performOperationWithName('replacement');
            testCase.verifyEqual([sharedA sharedB],[12 14])
            testCase.verifyEqual(replacement.callCount,1)
            testCase.verifyEqual(unrelated.callCount,1)

            early = WVCountingOperation('batch early',[TestOperationRegistrationAndCaching.annotation("batch_shared") TestOperationRegistrationAndCaching.annotation("early_only")],@(~) {20,21});
            late = WVCountingOperation('batch late',[TestOperationRegistrationAndCaching.annotation("batch_shared") TestOperationRegistrationAndCaching.annotation("late_only")],@(~) {30,31});
            wvt.addOperation([early late],shouldOverwriteExisting=true,shouldSuppressWarning=true);
            testCase.verifyError(@()wvt.operationWithName('batch early'),'')
            testCase.verifyTrue(wvt.operationWithName('batch late') == late)
            testCase.verifyEqual(wvt.batch_shared,30)
        end

        function removalUsesIdentityAndClearsOnlyOwnedState(testCase)
            wvt = TestOperationRegistrationAndCaching.transform();
            removableAnnotations = [TestOperationRegistrationAndCaching.annotation("remove_a") TestOperationRegistrationAndCaching.annotation("remove_b")];
            removable = WVCountingOperation('removable',removableAnnotations,@(~) {41,42});
            retained = WVCountingOperation('retained',TestOperationRegistrationAndCaching.annotation("retained"),@(~) {43});
            wvt.addOperation([removable retained]);
            wvt.performOperationWithName('removable');
            testCase.verifyEqual(wvt.retained,43)

            foreignAnnotations = [TestOperationRegistrationAndCaching.annotation("remove_a") TestOperationRegistrationAndCaching.annotation("remove_b")];
            foreign = WVCountingOperation('removable',foreignAnnotations,@(~) {51,52});
            testCase.verifyError(@()wvt.removeOperation(foreign),'')
            testCase.verifyTrue(wvt.operationWithName('removable') == removable)
            testCase.verifyTrue(all(isKey(wvt.variableCache,["remove_a" "remove_b" "retained"])))

            wvt.removeOperation(removable);
            testCase.verifyError(@()wvt.removeOperation(removable),'')
            testCase.verifyError(@()wvt.operationWithName('removable'),'')
            testCase.verifyError(@()wvt.performOperationWithName('removable'),'')
            testCase.verifyError(@()wvt.variableWithName('remove_a'),'')
            testCase.verifyFalse(any(ismember(["remove_a" "remove_b"],wvt.variableNames)))
            testCase.verifyFalse(any(isKey(wvt.variableCache,["remove_a" "remove_b"])))
            testCase.verifyTrue(isKey(wvt.variableCache,"retained"))
            testCase.verifyFalse(any(isKey(wvt.timeDependentVariablesNameMap,["remove_a" "remove_b"])))
            testCase.verifyFalse(any(isKey(wvt.wvCoefficientDependentVariablesNameMap,["remove_a" "remove_b"])))
        end

        function dependencyMetadataInvalidatesOnlyDeclaredOutputs(testCase)
            wvt = TestOperationRegistrationAndCaching.transform();
            timeAnnotation = TestOperationRegistrationAndCaching.annotation("time_only");
            timeAnnotation.isVariableWithLinearTimeStep = true;
            timeAnnotation.isDependentOnApAmA0 = false;
            coefficientAnnotation = TestOperationRegistrationAndCaching.annotation("coefficient_only");
            coefficientAnnotation.isVariableWithLinearTimeStep = false;
            coefficientAnnotation.isDependentOnApAmA0 = true;
            pair = WVCountingOperation('dependency pair',[timeAnnotation coefficientAnnotation],@(transform) {transform.t,sum(transform.Ap(:))+sum(transform.Am(:))+sum(transform.A0(:))});

            invariantAnnotation = TestOperationRegistrationAndCaching.annotation("invariant");
            invariantAnnotation.isVariableWithLinearTimeStep = false;
            invariantAnnotation.isDependentOnApAmA0 = false;
            invariant = WVCountingOperation('invariant',invariantAnnotation,@(~) {61});
            wvt.addOperation([pair invariant]);
            wvt.performOperationWithName('dependency pair');
            testCase.verifyEqual(wvt.invariant,61)
            testCase.verifyEqual(pair.callCount,1)
            testCase.verifyTrue(isKey(wvt.timeDependentVariablesNameMap,"time_only"))
            testCase.verifyFalse(isKey(wvt.timeDependentVariablesNameMap,"coefficient_only"))
            testCase.verifyTrue(isKey(wvt.wvCoefficientDependentVariablesNameMap,"coefficient_only"))
            testCase.verifyFalse(isKey(wvt.wvCoefficientDependentVariablesNameMap,"time_only"))

            wvt.t = wvt.t+1;
            testCase.verifyFalse(isKey(wvt.variableCache,"time_only"))
            testCase.verifyTrue(all(isKey(wvt.variableCache,["coefficient_only" "invariant"])))
            testCase.verifyEqual(wvt.coefficient_only,0)
            testCase.verifyEqual(pair.callCount,1)
            testCase.verifyEqual(wvt.time_only,wvt.t)
            testCase.verifyEqual(pair.callCount,2)

            coefficientProperties = ["Ap" "Am" "A0"];
            for iProperty = 1:length(coefficientProperties)
                propertyName = coefficientProperties(iProperty);
                values = wvt.(propertyName);
                values(2) = values(2)+iProperty;
                wvt.(propertyName) = values;
                testCase.verifyFalse(isKey(wvt.variableCache,"coefficient_only"))
                testCase.verifyTrue(all(isKey(wvt.variableCache,["time_only" "invariant"])))
                previousCount = pair.callCount;
                wvt.coefficient_only;
                testCase.verifyEqual(pair.callCount,previousCount+1)
            end

            replacementAnnotation = TestOperationRegistrationAndCaching.annotation("coefficient_only");
            replacementAnnotation.isVariableWithLinearTimeStep = true;
            replacementAnnotation.isDependentOnApAmA0 = false;
            replacement = WVCountingOperation('coefficient_only',replacementAnnotation,@(transform) {transform.t+100});
            wvt.addOperation(replacement,shouldOverwriteExisting=true,shouldSuppressWarning=true);
            testCase.verifyFalse(isKey(wvt.wvCoefficientDependentVariablesNameMap,"coefficient_only"))
            testCase.verifyTrue(isKey(wvt.timeDependentVariablesNameMap,"coefficient_only"))
            testCase.verifyEqual(wvt.coefficient_only,wvt.t+100)
        end
    end

    methods (Static, Access=private)
        function wvt = transform()
            wvt = WVTransformConstantStratification([4e3 3e3 2e3],[8 6 5],latitude=30,shouldAntialias=false);
        end

        function annotation = annotation(name)
            annotation = WVVariableAnnotation(char(name),{},'1','test operation output');
        end

        function state = registryState(wvt)
            state.operationNames = sort(string(wvt.operationNameMap.keys));
            state.variableNames = sort(string(wvt.operationVariableNameMap.keys));
            state.annotationNames = sort(string(wvt.annotatedPropertyNames));
            state.cacheNames = sort(string(wvt.variableCache.keys));
            state.coefficientDependentNames = sort(string(wvt.wvCoefficientDependentVariablesNameMap.keys));
            state.timeDependentNames = sort(string(wvt.timeDependentVariablesNameMap.keys));
        end

        function verifyRegistryState(testCase,wvt,expected)
            actual = TestOperationRegistrationAndCaching.registryState(wvt);
            testCase.verifyEqual(actual,expected)
        end
    end
end
