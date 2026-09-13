classdef WVTotalFlowComponent < WVPrimaryFlowComponent
    %Orthogonal solution group
    %
    % Each degree-of-freedom in the model is associated with an analytical
    % solution to the equations of motion. This class groups together
    % solutions of a particular type and provides a mapping between their
    % analytical solutions and their numerical representation.
    %
    % Perhaps the most complicate part of the numerical implementation is
    % the indexing---finding where each solution is represented
    % numerically. In general, a solution will have some properties, e.g.,
    %   (kMode,lMode,jMode,phi,A,omegasign) 
    % which will have a primary and conjugate part, each of which might be
    % in two different matrices.
    %
    % - Topic: Initialization
    % - Topic: Properties
    % - Topic: Masks
    % - Topic: Quadratic quantities
    % - Topic: Valid mode indices
    % - Topic: Analytical solutions
    properties (Access=private)
        bitmask = 0
    end

    methods
        function self = WVTotalFlowComponent(wvt)
            % create a new orthogonal solution group
            %
            % - Topic: Initialization
            % - Declaration:  solnGroup = WVFlowComponent(wvt)
            % - Parameter wvt: instance of a WVTransform
            % - Returns solnGroup: a new orthogonal solution group instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVPrimaryFlowComponent(wvt);
            self.name = "total";
            self.abbreviatedName = "tot";
        end

        function mask = maskOfPrimaryModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Masks
            % - Declaration: mask = maskOfPrimaryModesForCoefficientMatrix(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = 0;
            for i=1:length(self.wvt.primaryFlowComponents)
                mask = mask | self.wvt.primaryFlowComponents(i).maskOfPrimaryModesForCoefficientMatrix(coefficientMatrix);
            end
        end

        function mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Masks
            % - Declaration: mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = 0;
            for i=1:length(self.wvt.primaryFlowComponents)
                mask = mask | self.wvt.primaryFlowComponents(i).maskOfConjugateModesForCoefficientMatrix(coefficientMatrix);
            end
        end

        function totalEnergyFactor = totalEnergyFactorForCoefficientMatrix(self,coefficientMatrix)
            % returns the total energy multiplier for the coefficient matrix.
            %
            % Returns a matrix of size wvt.spectralMatrixSize that
            % multiplies the squared absolute value of this matrix to
            % produce the total energy.
            %
            % - Topic: Quadratic quantities
            % - Declaration: totalEnergyFactor = totalEnergyFactorForCoefficientMatrix(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                totalEnergyFactor double {mustBeNonnegative}
            end
            totalEnergyFactor = zeros(self.wvt.spectralMatrixSize);
            for i=1:length(self.wvt.primaryFlowComponents)
                totalEnergyFactor = totalEnergyFactor + self.wvt.primaryFlowComponents(i).totalEnergyFactor(coefficientMatrix);
            end
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)
            % returns a boolean indicating whether (k,l,j) is a valid mode number
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % number
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns index: a non-negative integer
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            for i=1:length(self.wvt.primaryFlowComponents)
                bool = bool | self.wvt.primaryFlowComponents(i).isValidPrimaryModeNumber(kMode,lMode,jMode);
            end
        end

        function bool = isValidConjugateModeNumber(self,kMode,lMode,jMode)
            % returns a boolean indicating whether (k,l,j) is a valid mode number
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % number
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns index: a non-negative integer
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            for i=1:length(self.wvt.primaryFlowComponents)
                bool = bool | self.wvt.primaryFlowComponents(i).isValidConjugateModeNumber(kMode,lMode,jMode);
            end
        end

        function solutions = solutionForModeAtIndex(self,index,options)
            % Return analytical solutions from the complete primary-flow basis.
            %
            % Total-flow indices run from 1 through `nModes`. Primary flow
            % components are ordered lexically by `shortName`, and each
            % component contributes one contiguous range while retaining its
            % own local mode ordering. Scalar inputs return a scalar solution;
            % column-vector inputs return solutions in the requested order.
            %
            % Set `amplitude='wvt'` to reconstruct each solution from the
            % corresponding coefficient currently stored by the transform.
            % Set `amplitude='random'` to generate an appropriate random
            % amplitude for each requested solution.
            %
            % - Topic: Analytical solutions
            % - Declaration: solutions = solutionForModeAtIndex(index,options)
            % - Parameter index: positive integer scalar or column vector with values no greater than nModes
            % - Parameter amplitude: (optional) 'wvt' or 'random' (default)
            % - Returns solutions: scalar or column vector of WVOrthogonalSolution objects
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                index (:,1) double {mustBeInteger,mustBePositive}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            mustBeLessThanOrEqual(index,self.nModes)

            components = self.wvt.primaryFlowComponents;
            [~,componentOrder] = sort(string({components.shortName}));
            components = components(componentOrder);
            nModesByComponent = arrayfun(@(component) component.nModes,components);
            lastIndexByComponent = cumsum(nModesByComponent);

            solutions = WVOrthogonalSolution.empty(length(index),0);
            for iSolution = 1:length(index)
                iComponent = find(index(iSolution) <= lastIndexByComponent,1);
                firstIndexForComponent = lastIndexByComponent(iComponent) - nModesByComponent(iComponent) + 1;
                localIndex = index(iSolution) - firstIndexForComponent + 1;
                solutions(iSolution) = components(iComponent).solutionForModeAtIndex(localIndex,amplitude=options.amplitude);
            end
        end


    end
    
end
