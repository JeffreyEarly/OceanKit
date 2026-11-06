classdef (Hidden) IMHydrostaticInnerProductCatalog
    % Map hydrostatic boundary conditions to endpoint inner-product terms.

    methods (Static)
        function result = resolve(evp, variable)
            arguments
                evp IMInternalModes
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])}
            end

            variable = string(variable);
            result.hasInnerProduct = false;
            result.reason = "The diagnostic inner product is available only for " ...
                + "modeFamily=""hydrostatic"" EVPs in the value-only catalog.";
            result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();

            if evp.modeFamily ~= "hydrostatic"
                return;
            end
            if variable == evp.formulation
                result.hasInnerProduct = true;
                result.reason = "The solved formulation uses the canonical scalar EVP inner product.";
                return;
            end

            surface = IMHydrostaticInnerProductCatalog.endpointResult(evp, variable, "surface", evp.surfaceBoundary);
            bottom = IMHydrostaticInnerProductCatalog.endpointResult(evp, variable, "bottom", evp.bottomBoundary);
            result.endpointInnerProductTerms = [surface.endpointInnerProductTerms; bottom.endpointInnerProductTerms];
            result.hasInnerProduct = surface.hasInnerProduct && bottom.hasInnerProduct;
            if result.hasInnerProduct && isempty(result.endpointInnerProductTerms)
                result.reason = "The hydrostatic value-only catalog gives an " ...
                    + "interior-only diagnostic inner product for this boundary combination.";
            elseif result.hasInnerProduct
                result.reason = "The hydrostatic value-only catalog gives a " ...
                    + "known diagnostic inner product with endpoint inner-product terms.";
            else
                result.reason = IMHydrostaticInnerProductCatalog.joinReasons([surface.reason bottom.reason]);
            end
        end

        function terms = emptyEndpointInnerProductTerms()
            terms = struct("location", {}, "coefficient", {}, ...
                "variable", {}, "description", {}, "catalogCase", {});
        end
    end

    methods (Static, Access = private)
        function result = endpointResult(evp, variable, location, boundary)
            coefficients = [boundary.a boundary.b boundary.c boundary.d];
            tolerance = 100*eps*max([1 abs(coefficients)]);
            active = abs(coefficients) > tolerance;
            result.hasInnerProduct = true;
            result.reason = "This endpoint adds no diagnostic endpoint term.";
            result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();

            if nnz(active) <= 1
                return;
            end
            if nnz(active) == 4
                result.hasInnerProduct = false;
                result.reason = "This boundary condition is outside the value-only diagnostic inner-product catalog.";
                return;
            end

            if evp.formulation == "F" && variable == "G"
                result = IMHydrostaticInnerProductCatalog.fFormulationGEndpoint( ...
                    location, boundary, active, evp.g);
            elseif evp.formulation == "G" && variable == "F"
                result = IMHydrostaticInnerProductCatalog.gFormulationFEndpoint(location, boundary, active);
            else
                result.hasInnerProduct = false;
                result.reason = "No diagnostic inner-product catalog entry is " ...
                    + "available for this formulation and variable.";
            end
        end

        function result = fFormulationGEndpoint(location, boundary, active, g)
            a = boundary.a;
            b = boundary.b;
            c = boundary.c;
            d = boundary.d;
            eta = IMHydrostaticInnerProductCatalog.eta(location);
            result.hasInnerProduct = true;
            result.reason = "The hydrostatic catalog provides a value-only diagnostic G inner product.";
            result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();

            if isequal(active, [false false true true])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*d/c, "G", "F-P1");
            elseif isequal(active, [true true false false])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*b/(g*a), "G", "F-P2");
            elseif isequal(active, [false true true true])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*d/c, "G", "F-T1");
            elseif isequal(active, [true true true false])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*g*a/b, "F", "F-T4");
            elseif isequal(active, [false true false true]) || isequal(active, [true false true false]) ...
                    || isequal(active, [false true true false])
                result.reason = "The diagnostic G inner product is the interior integral for this endpoint.";
            elseif isequal(active, [true false false true]) || isequal(active, [true true false true])
                result.hasInnerProduct = false;
                result.reason = "This boundary condition gives only a mixed identity, " ...
                    + "not a standalone diagnostic inner product.";
            elseif isequal(active, [true false true true])
                result.hasInnerProduct = false;
                result.reason = "This boundary condition gives eigenvalue-dependent endpoint terms, " ...
                    + "which are not supported as a known diagnostic inner product.";
            elseif isequal(active, [true false false false]) || isequal(active, [false true false false]) ...
                    || isequal(active, [false false true false]) || isequal(active, [false false false true])
                result.reason = "This endpoint adds no diagnostic endpoint term.";
            else
                result.hasInnerProduct = false;
                result.reason = "This boundary condition does not provide a " ...
                    + "known value-only diagnostic inner product.";
            end
        end

        function result = gFormulationFEndpoint(location, boundary, active)
            a = boundary.a;
            b = boundary.b;
            c = boundary.c;
            d = boundary.d;
            eta = IMHydrostaticInnerProductCatalog.eta(location);
            result.hasInnerProduct = true;
            result.reason = "The hydrostatic catalog provides a value-only diagnostic F inner product.";
            result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();

            if isequal(active, [true true false false])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*b/a, "F", "G-P5");
            elseif isequal(active, [false false true true])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*d/c, "F", "G-P6");
            elseif isequal(active, [false true true true])
                result.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.innerProductTerm( ...
                    location, -eta*d/c, "F", "G-T4");
            elseif isequal(active, [false true true false]) || isequal(active, [true false true false]) ...
                    || isequal(active, [false true false true])
                result.reason = "The diagnostic F inner product is the interior integral for this endpoint.";
            elseif isequal(active, [true false false true]) || isequal(active, [true false true true]) ...
                    || isequal(active, [true true false true])
                result.hasInnerProduct = false;
                result.reason = "This boundary condition gives eigenvalue-dependent endpoint terms, " ...
                    + "which are not supported as a known diagnostic inner product.";
            elseif isequal(active, [true true true false])
                result.hasInnerProduct = false;
                result.reason = "This boundary condition requires derivative endpoint terms, which are not supported.";
            elseif isequal(active, [true false false false]) || isequal(active, [false false true false]) ...
                    || isequal(active, [false false false true])
                result.reason = "This endpoint adds no diagnostic endpoint term.";
            else
                result.hasInnerProduct = false;
                result.reason = "This boundary condition does not provide a " ...
                    + "known value-only diagnostic inner product.";
            end
        end

        function term = innerProductTerm(location, coefficient, variable, catalogCase)
            term = struct( ...
                "location", string(location), ...
                "coefficient", coefficient, ...
                "variable", string(variable), ...
                "description", string(variable) + "(" + string(location) + ")", ...
                "catalogCase", string(catalogCase));
        end

        function eta = eta(location)
            if string(location) == "surface"
                eta = 1;
            else
                eta = -1;
            end
        end

        function reason = joinReasons(reasons)
            reasons = unique(string(reasons), "stable");
            reasons(strlength(reasons) == 0) = [];
            if isempty(reasons)
                reason = "";
            else
                reason = join(reasons, " ");
            end
        end
    end
end
