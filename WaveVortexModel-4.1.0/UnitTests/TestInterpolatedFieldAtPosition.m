classdef TestInterpolatedFieldAtPosition < matlab.unittest.TestCase

    properties
        wvt
        field2D
        field3D
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            Lxyz = [1000 800 500];
            Nxyz = [8 8 7];
            N2 = @(z) (5.2e-3)^2*(1 + 0.2*z/Lxyz(3));
            testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,shouldAntialias=false);

            X = testCase.wvt.X;
            Y = testCase.wvt.Y;
            Z = testCase.wvt.Z;
            surfaceField = 0.7*sin(2*pi*X(:,:,end)/testCase.wvt.Lx) ...
                - 0.3*cos(2*pi*Y(:,:,end)/testCase.wvt.Ly) ...
                + 0.2*sin(2*pi*(X(:,:,end)/testCase.wvt.Lx + Y(:,:,end)/testCase.wvt.Ly));
            volumeField = surfaceField.*(1 + 0.2*(Z/testCase.wvt.Lz).^2) + 1e-3*Z;

            testCase.field2D = surfaceField;
            testCase.field3D = volumeField;

            annotation2D = WVVariableAnnotation('interpolationTest2D',{'x','y'},'1','synthetic two-dimensional interpolation test field');
            annotation3D = WVVariableAnnotation('interpolationTest3D',{'x','y','z'},'1','synthetic three-dimensional interpolation test field');
            testCase.wvt.addOperation(WVOperation('interpolationTest2D',annotation2D,@(~) surfaceField));
            testCase.wvt.addOperation(WVOperation('interpolationTest3D',annotation3D,@(~) volumeField));
        end
    end

    methods (Test)
        function testTwoDimensionalLinearInterpolation(testCase)
            [x,y] = testCase.horizontalQueryPoints();

            actual = testCase.wvt.variableAtPositionWithName(x,y,[],'interpolationTest2D');
            expected = testCase.periodicInterpolation2D(testCase.field2D,x,y);

            testCase.verifySize(actual,size(x));
            testCase.verifyEqual(actual,expected,AbsTol=1e-12);
        end

        function testThreeDimensionalLinearInterpolation(testCase)
            [x,y] = testCase.horizontalQueryPoints();
            zGrid = testCase.wvt.z;
            testCase.verifyGreaterThan(max(diff(zGrid))-min(diff(zGrid)),1e-6);
            z = reshape((zGrid(2:5)+zGrid(3:6))/2,1,[]);

            actual = testCase.wvt.variableAtPositionWithName(x,y,z,'interpolationTest3D');
            expected = testCase.periodicInterpolation3D(testCase.field3D,x,y,z);

            testCase.verifySize(actual,size(x));
            testCase.verifyEqual(actual,expected,AbsTol=1e-12);
        end
    end

    methods (Access=private)
        function [x,y] = horizontalQueryPoints(testCase)
            dx = testCase.wvt.x(2)-testCase.wvt.x(1);
            dy = testCase.wvt.y(2)-testCase.wvt.y(1);

            % interior, x boundary, y boundary, x-y corner
            x = [3.25*dx -0.25*dx 3.25*dx testCase.wvt.Lx+0.2*dx];
            y = [3.4*dy 3.4*dy testCase.wvt.Ly+0.3*dy testCase.wvt.Ly-0.25*dy];
        end

        function expected = periodicInterpolation2D(testCase,field,x,y)
            [xPeriodic,yPeriodic] = testCase.periodicGrid();
            fieldPeriodic = field([end 1:end 1],[end 1:end 1]);
            expected = interpn(xPeriodic,yPeriodic,fieldPeriodic,mod(x,testCase.wvt.Lx),mod(y,testCase.wvt.Ly),'linear');
        end

        function expected = periodicInterpolation3D(testCase,field,x,y,z)
            [xPeriodic,yPeriodic] = testCase.periodicGrid();
            fieldPeriodic = field([end 1:end 1],[end 1:end 1],:);
            expected = interpn(xPeriodic,yPeriodic,testCase.wvt.z,fieldPeriodic,mod(x,testCase.wvt.Lx),mod(y,testCase.wvt.Ly),z,'linear');
        end

        function [xPeriodic,yPeriodic] = periodicGrid(testCase)
            xPeriodic = [testCase.wvt.x(end)-testCase.wvt.Lx;testCase.wvt.x;testCase.wvt.x(1)+testCase.wvt.Lx];
            yPeriodic = [testCase.wvt.y(end)-testCase.wvt.Ly;testCase.wvt.y;testCase.wvt.y(1)+testCase.wvt.Ly];
        end
    end

end
