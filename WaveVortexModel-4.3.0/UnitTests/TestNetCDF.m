classdef TestNetCDF < matlab.unittest.TestCase
    properties
        path
        x = 0:9
        f_a = @(x) 4*x + sqrt(-1)*2*x
        b = logical([0 1 0 1 1 0 1 0 1 1])
    end

    methods (TestMethodSetup)
        function createTemporaryPath(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.path = fullfile(fixture.Folder,'test.nc');
        end
    end

    methods (Test, TestTags = "full")
        function testAddDimension(testCase)
            ncfile = testCase.createFile();
            testCase.verifyWarningFree(@()ncfile.addDimension('x',testCase.x))
        end

        function testReadDimension(testCase)
            ncfile = testCase.createFile();
            ncfile.addDimension('x',testCase.x);
            xBack = ncfile.readVariables('x');
            testCase.verifyEqual(xBack,testCase.x)
        end

        function testAddComplexVariable(testCase)
            ncfile = testCase.createFileWithXDimension();
            testCase.verifyWarningFree(@()ncfile.addVariable('a',{'x'},testCase.f_a(testCase.x)))
        end

        function testReadComplexVariable(testCase)
            ncfile = testCase.createFileWithXDimension();
            ncfile.addVariable('a',{'x'},testCase.f_a(testCase.x));
            aBack = ncfile.readVariables('a');
            testCase.verifyEqual(aBack,testCase.f_a(testCase.x))
        end

        function testAddLogicalVariable(testCase)
            ncfile = testCase.createFileWithXDimension();
            testCase.verifyWarningFree(@()ncfile.addVariable('b',{'x'},testCase.b))
        end

        function testReadLogicalVariable(testCase)
            ncfile = testCase.createFileWithXDimension();
            ncfile.addVariable('b',{'x'},testCase.b);
            bBack = ncfile.readVariables('b');
            testCase.verifyEqual(bBack,testCase.b)
        end

        function testAddColumnVector(testCase)
            ncfile = testCase.createFileWithXDimension();
            testCase.verifyWarningFree(@()ncfile.addVariable('y',{'x'},reshape(testCase.x,[],1)))
        end

        function testReadColumnVector(testCase)
            ncfile = testCase.createFileWithXDimension();
            y = reshape(testCase.x,[],1);
            ncfile.addVariable('y',{'x'},y);
            yBack = ncfile.readVariables('y');
            testCase.verifyEqual(yBack,y)
            testCase.verifyNotEqual(yBack,testCase.x)
        end

        function testUnlimitedDimension(testCase)
            ncfile = testCase.createFile();
            s = [1 2.5 3.14].';
            [~,variable] = ncfile.addDimension('t',length=Inf,type='double');
            for i=1:length(s)
                variable.setValueAlongDimensionAtIndex(s(i),{'t'},i);
            end

            for i=1:length(s)
                testCase.verifyEqual(variable.valueAlongDimensionAtIndex({'t'},i),s(i))
                testCase.verifyEqual(ncfile.readVariablesAtIndexAlongDimension({'t'},i,'t'),s(i))
            end
            testCase.verifyEqual(variable.value,s)
        end

        function testAddGroup(testCase)
            ncfile = testCase.createFile();
            group = ncfile.addGroup("MyGroup");
            group.addAttribute('MyAttribute',"Hello group!")
            s = 0:9;
            group.addDimension('s',s);
            values = 4*s + sqrt(-1)*2*s;
            group.addVariable('AFunction',{'s'},values);

            valuesBack = ncfile.readVariables('MyGroup/AFunction');
            testCase.verifyEqual(valuesBack,values)
        end
    end

    methods (Access=private)
        function ncfile = createFile(testCase)
            ncfile = NetCDFFile(testCase.path,shouldOverwriteExisting=true);
            testCase.addTeardown(@()TestNetCDF.closeIfOpen(ncfile));
        end

        function ncfile = createFileWithXDimension(testCase)
            ncfile = testCase.createFile();
            ncfile.addDimension('x',testCase.x);
        end
    end

    methods (Static, Access=private)
        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
