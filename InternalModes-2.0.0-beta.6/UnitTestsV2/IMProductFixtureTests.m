classdef IMProductFixtureTests < matlab.unittest.TestCase
    methods (Test)
        function frozenMixedProductControlsAgree(testCase)
            folder=fullfile(fileparts(mfilename("fullpath")),"Fixtures");
            fixture=jsondecode(fileread(fullfile(folder,"wvm-cal-constant-17-products.json")));
            grids=reshape(fixture.grids,1,[]);
            for j=1:numel(grids), grids(j).id=string(grids(j).id); end
            entries=cellRow(fixture.factors); factors=cell(size(entries));
            for j=1:numel(entries), factors{j}=factorFromFixture(entries{j},grids); end
            entries=cellRow(fixture.outputs); outputs=cell(size(entries));
            for j=1:numel(entries), outputs{j}=outputFromFixture(entries{j}); end
            products=struct2table(fixture.products);
            result=IMProductAssessmentPlan(factors,products,outputs,retainedCounts=fixture.prefixCounts.',productBudget=6).assess(grids,chunkSize=2,minimumReciprocalCondition=1e-13);
            testCase.verifyClass(result,"struct")
            testCase.verifyEqual(result.measurements.retainedCount,fixture.prefixCounts(:))
            errors=cellRow(fixture.expectedErrors);
            for j=1:numel(errors)
                testCase.verifyEqual(result.evidence{j}.error,unpack(errors{j}),AbsTol=2e-12)
                testCase.verifyEqual(result.evidence{j}.isZero,logical(fixture.expectedZeros(j)))
                testCase.verifyEqual(result.evidence{j}.referenceStatus,"qualified")
            end
            testCase.verifyEqual(result.costs.reservedProducts,6)
            testCase.verifyEqual(unique(string(cellfun(@(f) f.family,factors,UniformOutput=false))),["apv","boundary","wave"])
            testCase.verifyEqual(fixture.verification.comparisonCount,639744)
            testCase.verifyLessThan(fixture.verification.maximumSinglePrecisionErrorDifference,1e-12)
        end
    end
end

function factor=factorFromFixture(factor,grids)
samples=cellfun(@unpack,cellRow(factor.samples),UniformOutput=false);
factor=rmfield(factor,"samples");
factor.evaluate=@evaluate;
    function values=evaluate(z,columns,id)
        index=find([grids.id]==id);
        assert(isscalar(index) && isequal(z,grids(index).z),'Fixture grid identity must be preserved.');
        values=samples{index}(:,columns);
    end
end

function output=outputFromFixture(output)
c=output.context; output=rmfield(output,"context");
projection=IMProjection.fromPrescribedDual(unpack(c.samplePairingMatrix),c.sampleGram,c.targetGram,majorantGramMatrix=c.majorantGram,activeColumnMask=reshape(c.active,1,[]),columnLabels=reshape(string(c.labels),1,[]),provenance=c.provenance);
references=cellRow(c.references);
for j=1:numel(references)
    reference=references{j};
    reference.pairingMatrix=unpack(reference.pairingMatrix);
    reference.normMatrix=diag(reference.normDiagonal);
    references{j}=rmfield(reference,"normDiagonal");
end
context=struct(projection=projection,references={references});
output.prepare=@(grids) context;
end

function cells=cellRow(values)
if isstruct(values), cells=num2cell(values); else, cells=values; end
cells=reshape(cells,1,[]);
end

function values=unpack(packed)
realValues=zeros(prod(packed.shape),1); imaginaryValues=realValues;
if ~isempty(packed.real), realValues(:)=packed.real; end
if ~isempty(packed.imag), imaginaryValues(:)=packed.imag; end
values=reshape(complex(realValues,imaginaryValues),packed.shape.');
end
