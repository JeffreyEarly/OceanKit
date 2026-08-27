classdef TestInterpolatedFieldAtPosition < matlab.unittest.TestCase

    properties
        wvt
        barotropicTransform
        field2D
        field3D
        barotropicField
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

            testCase.barotropicTransform = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=30,shouldAntialias=false);
            XBarotropic = testCase.barotropicTransform.X;
            YBarotropic = testCase.barotropicTransform.Y;
            testCase.barotropicField = 0.4*cos(2*pi*XBarotropic/testCase.barotropicTransform.Lx) ...
                + 0.6*sin(2*pi*YBarotropic/testCase.barotropicTransform.Ly) ...
                - 0.1*cos(2*pi*(XBarotropic/testCase.barotropicTransform.Lx-YBarotropic/testCase.barotropicTransform.Ly));
            barotropicAnnotation = WVVariableAnnotation('barotropicInterpolationTest',{'x','y'},'1','synthetic barotropic interpolation test field');
            testCase.barotropicTransform.addOperation(WVOperation('barotropicInterpolationTest',barotropicAnnotation,@(~) testCase.barotropicField));
        end
    end

    methods (Test, TestTags = "smoke")
        function testTwoDimensionalLinearInterpolation(testCase)
            [x,y] = testCase.horizontalQueryPoints();

            actual = testCase.wvt.variableAtPositionWithName(x,y,[],'interpolationTest2D');
            expected = testCase.periodicInterpolation2D(testCase.wvt,testCase.field2D,x,y);

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

        function testSupportedMethodsReturnMultipleFields(testCase)
            [x,y] = testCase.horizontalQueryPoints();
            zGrid = testCase.wvt.z;
            z = reshape((zGrid(2:5)+zGrid(3:6))/2,1,[]);

            for method = ["linear" "spline"]
                [actual2D,actual3D] = testCase.wvt.variableAtPositionWithName(x,y,z,'interpolationTest2D','interpolationTest3D',interpolationMethod=method);
                expected2D = testCase.wvt.variableAtPositionWithName(x,y,z,'interpolationTest2D',interpolationMethod=method);
                expected3D = testCase.wvt.variableAtPositionWithName(x,y,z,'interpolationTest3D',interpolationMethod=method);

                testCase.verifySize(actual2D,size(x));
                testCase.verifySize(actual3D,size(x));
                testCase.verifyEqual(actual2D,expected2D);
                testCase.verifyEqual(actual3D,expected3D);
            end
        end

        function testSupportedMethodsArePeriodic(testCase)
            [x,y] = testCase.horizontalQueryPoints();
            zGrid = testCase.wvt.z;
            z = reshape((zGrid(2:5)+zGrid(3:6))/2,1,[]);
            xShifted = x + [testCase.wvt.Lx -testCase.wvt.Lx 2*testCase.wvt.Lx -2*testCase.wvt.Lx];
            yShifted = y + [-testCase.wvt.Ly 2*testCase.wvt.Ly -2*testCase.wvt.Ly testCase.wvt.Ly];

            for method = ["linear" "spline"]
                [expected2D,expected3D] = testCase.wvt.variableAtPositionWithName(x,y,z,'interpolationTest2D','interpolationTest3D',interpolationMethod=method);
                [actual2D,actual3D] = testCase.wvt.variableAtPositionWithName(xShifted,yShifted,z,'interpolationTest2D','interpolationTest3D',interpolationMethod=method);

                testCase.verifySize(actual2D,size(x));
                testCase.verifySize(actual3D,size(x));
                testCase.verifyEqual(actual2D,expected2D,AbsTol=1e-12);
                testCase.verifyEqual(actual3D,expected3D,AbsTol=1e-12);
            end
        end

        function testSupportedMethodsReproduceGridValues(testCase)
            iX = [1 4 8];
            iY = [2 5 1];
            iZ = [1 4 7];
            x = reshape(testCase.wvt.x(iX),1,[]);
            y = reshape(testCase.wvt.y(iY),1,[]);
            z = reshape(testCase.wvt.z(iZ),1,[]);
            expected2D = reshape(testCase.field2D(sub2ind(size(testCase.field2D),iX,iY)),1,[]);
            expected3D = reshape(testCase.field3D(sub2ind(size(testCase.field3D),iX,iY,iZ)),1,[]);

            for method = ["linear" "spline"]
                [actual2D,actual3D] = testCase.wvt.variableAtPositionWithName(x,y,z,'interpolationTest2D','interpolationTest3D',interpolationMethod=method);
                testCase.verifyEqual(actual2D,expected2D,AbsTol=1e-12);
                testCase.verifyEqual(actual3D,expected3D,AbsTol=1e-12);
            end
        end

        function testSplineIsContinuousAcrossPeriodicSeams(testCase)
            dx = testCase.wvt.x(2)-testCase.wvt.x(1);
            dy = testCase.wvt.y(2)-testCase.wvt.y(1);
            epsilonX = 1e-7*dx;
            epsilonY = 1e-7*dy;
            z = testCase.wvt.z(4);

            [left2D,left3D] = testCase.wvt.variableAtPositionWithName(testCase.wvt.Lx-epsilonX,2.4*dy,z,'interpolationTest2D','interpolationTest3D',interpolationMethod='spline');
            [right2D,right3D] = testCase.wvt.variableAtPositionWithName(epsilonX,2.4*dy,z,'interpolationTest2D','interpolationTest3D',interpolationMethod='spline');
            [bottom2D,bottom3D] = testCase.wvt.variableAtPositionWithName(3.2*dx,testCase.wvt.Ly-epsilonY,z,'interpolationTest2D','interpolationTest3D',interpolationMethod='spline');
            [top2D,top3D] = testCase.wvt.variableAtPositionWithName(3.2*dx,epsilonY,z,'interpolationTest2D','interpolationTest3D',interpolationMethod='spline');

            tolerance = 1e-5*max(abs([testCase.field2D(:);testCase.field3D(:)]));
            testCase.verifyLessThan(abs(right2D-left2D),tolerance);
            testCase.verifyLessThan(abs(right3D-left3D),tolerance);
            testCase.verifyLessThan(abs(top2D-bottom2D),tolerance);
            testCase.verifyLessThan(abs(top3D-bottom3D),tolerance);
        end

        function testBarotropicSupportedInterpolation(testCase)
            transform = testCase.barotropicTransform;
            dx = transform.x(2)-transform.x(1);
            dy = transform.y(2)-transform.y(1);
            x = [3.25*dx -0.25*dx 3.25*dx transform.Lx+0.2*dx];
            y = [3.4*dy 3.4*dy transform.Ly+0.3*dy transform.Ly-0.25*dy];
            xShifted = x + [transform.Lx -transform.Lx 2*transform.Lx -2*transform.Lx];
            yShifted = y + [-transform.Ly 2*transform.Ly -2*transform.Ly transform.Ly];

            for method = ["linear" "spline"]
                expected = transform.variableAtPositionWithName(x,y,[],'barotropicInterpolationTest',interpolationMethod=method);
                actual = transform.variableAtPositionWithName(xShifted,yShifted,[],'barotropicInterpolationTest',interpolationMethod=method);
                testCase.verifySize(actual,size(x));
                testCase.verifyEqual(actual,expected,AbsTol=1e-12);
            end

            expectedLinear = testCase.periodicInterpolation2D(transform,testCase.barotropicField,x,y);
            actualLinear = transform.variableAtPositionWithName(x,y,[],'barotropicInterpolationTest',interpolationMethod='linear');
            testCase.verifyEqual(actualLinear,expectedLinear,AbsTol=1e-12);

            iX = [1 4 8];
            iY = [2 5 1];
            xGrid = reshape(transform.x(iX),1,[]);
            yGrid = reshape(transform.y(iY),1,[]);
            expectedGrid = reshape(testCase.barotropicField(sub2ind(size(testCase.barotropicField),iX,iY)),1,[]);
            for method = ["linear" "spline"]
                actualGrid = transform.variableAtPositionWithName(xGrid,yGrid,[],'barotropicInterpolationTest',interpolationMethod=method);
                testCase.verifyEqual(actualGrid,expectedGrid,AbsTol=1e-12);
            end

            epsilonX = 1e-7*dx;
            epsilonY = 1e-7*dy;
            left = transform.variableAtPositionWithName(transform.Lx-epsilonX,2.4*dy,[],'barotropicInterpolationTest',interpolationMethod='spline');
            right = transform.variableAtPositionWithName(epsilonX,2.4*dy,[],'barotropicInterpolationTest',interpolationMethod='spline');
            bottom = transform.variableAtPositionWithName(3.2*dx,transform.Ly-epsilonY,[],'barotropicInterpolationTest',interpolationMethod='spline');
            top = transform.variableAtPositionWithName(3.2*dx,epsilonY,[],'barotropicInterpolationTest',interpolationMethod='spline');
            tolerance = 1e-5*max(abs(testCase.barotropicField(:)));
            testCase.verifyLessThan(abs(right-left),tolerance);
            testCase.verifyLessThan(abs(top-bottom),tolerance);
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

        function expected = periodicInterpolation2D(~,transform,field,x,y)
            xPeriodic = [transform.x(end)-transform.Lx;transform.x;transform.x(1)+transform.Lx];
            yPeriodic = [transform.y(end)-transform.Ly;transform.y;transform.y(1)+transform.Ly];
            fieldPeriodic = field([end 1:end 1],[end 1:end 1]);
            expected = interpn(xPeriodic,yPeriodic,fieldPeriodic,mod(x,transform.Lx),mod(y,transform.Ly),'linear');
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
