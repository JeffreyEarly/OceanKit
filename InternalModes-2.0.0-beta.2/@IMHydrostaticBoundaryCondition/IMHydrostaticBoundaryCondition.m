classdef IMHydrostaticBoundaryCondition
    % Convert a hydrostatic `F`/`G` endpoint law to canonical form.
    %
    % `IMHydrostaticBoundaryCondition` stores one physical hydrostatic
    % endpoint law:
    %
    % $$
    % g\left[a+\frac{c}{gh}\right]F(z_\ell)
    % =
    % \left[b+\frac{d}{gh}+egh\right]G(z_\ell).
    % $$
    %
    % The coefficients `a`, `b`, `c`, `d`, and `e` belong to this
    % physical endpoint law. They are not the raw constructor coefficients
    % of `IMBoundaryCondition`.
    %
    % The known hydrostatic endpoint rules give endpoint additions to the
    % physical bilinear forms
    %
    % $$
    % \langle F_i,F_j\rangle_F
    % =
    % \int_{z_b}^{z_s}F_i(z)F_j(z)\,dz
    % +
    % \sum_\ell \Delta_F^\ell,
    % $$
    %
    % and
    %
    % $$
    % \langle G_i,G_j\rangle_G
    % =
    % \frac{1}{g}\int_{z_b}^{z_s}N^2(z)G_i(z)G_j(z)\,dz
    % +
    % \sum_\ell \Delta_G^\ell.
    % $$
    %
    % The table lists the endpoint additions $$\Delta_F^\ell$$ and
    % $$\Delta_G^\ell$$ for supported rows:
    %
    % | Endpoint law | $$\Delta_F^\ell$$ | $$\Delta_G^\ell$$ |
    % | --- | --- | --- |
    % | $$gaF=bG$$ | $$0$$ | $$\eta_\ell\frac{b}{ga}G^i_\ell G^j_\ell$$ |
    % | $$\frac{c}{h}F=bG$$ | $$-\eta_\ell\frac{c}{b}F^i_\ell F^j_\ell$$ | $$0$$ |
    % | $$g\left(a+\frac{c}{gh}\right)F=0$$ | $$0$$ | $$0$$ |
    % | $$g\left(a+\frac{c}{gh}\right)F=bG$$ | $$-\eta_\ell\frac{c}{b}F^i_\ell F^j_\ell$$ | $$\eta_\ell\frac{ga}{b}F^i_\ell F^j_\ell$$ |
    % | $$\frac{c}{h}F=\frac{d}{gh}G$$ | $$0$$ | $$\eta_\ell\frac{d}{gc}G^i_\ell G^j_\ell$$ |
    % | $$\left(b+\frac{d}{gh}\right)G=0$$ | $$0$$ | $$0$$ |
    % | $$(b+egh)G=0$$ | $$0$$ | $$0$$ |
    % | $$gaF=eghG$$ | $$-\eta_\ell\frac{a}{e}F^i_\ell F^j_\ell$$ | $$0$$ |
    % | $$\frac{c}{h}F=\left(b+\frac{d}{gh}\right)G$$ | $$-\eta_\ell\frac{1}{bc}\left(cF^i_\ell-\frac{d}{g}G^i_\ell\right)\left(cF^j_\ell-\frac{d}{g}G^j_\ell\right)$$ | $$\eta_\ell\frac{d}{gc}G^i_\ell G^j_\ell$$ |
    % | $$gaF=(b+egh)G$$ | $$-\eta_\ell\frac{1}{ae}\left(aF^i_\ell-\frac{b}{g}G^i_\ell\right)\left(aF^j_\ell-\frac{b}{g}G^j_\ell\right)$$ | $$\eta_\ell\frac{b}{ga}G^i_\ell G^j_\ell$$ |
    %
    % Conversion to canonical EVP boundaries is handled by
    % `canonicalBoundary`. For an `F`-formulated EVP, the law must have
    % $$e=0$$ and converts as
    %
    % ```matlab
    % boundary = IMBoundaryCondition(a=-law.a,b=law.b,c=-law.c/g,d=law.d/g);
    % ```
    %
    % For a `G`-formulated EVP, the law must have $$d=0$$ and converts as
    %
    % ```matlab
    % boundary = IMBoundaryCondition(a=law.e,b=law.a,c=law.b/g,d=law.c/g);
    % ```
    %
    % If the requested formulation would put both $$h$$ and $$1/h$$ on the
    % same side of the endpoint law, the conversion is nonlinear in
    % $$\lambda=1/h$$ and cannot be represented by one canonical linear
    % boundary row.
    %
    % After conversion, `IMInternalModes.innerProduct("F")` and
    % `IMInternalModes.innerProduct("G")` report which bilinear forms are
    % available for the resolved canonical boundary condition.
    %
    % ```matlab
    % g = 9.81;
    % law = IMHydrostaticBoundaryCondition(a=A/g,b=1);
    % surfaceBoundary = law.canonicalBoundary(formulation="G",g=g);
    % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0],g=g,surfaceBoundary=surfaceBoundary);
    % ```
    %
    % - Topic: Create hydrostatic boundary laws
    % - Topic: Convert hydrostatic boundary laws
    % - Topic: Inspect hydrostatic boundary laws
    % - Declaration: classdef IMHydrostaticBoundaryCondition

    properties (SetAccess = private)
        % Coefficient multiplying `F` on the \(h^0\) side.
        %
        % - Topic: Inspect hydrostatic boundary laws
        a = 0

        % Coefficient multiplying `G` on the \(h^0\) side.
        %
        % - Topic: Inspect hydrostatic boundary laws
        b = 0

        % Coefficient multiplying `F` on the \(1/h\) side.
        %
        % - Topic: Inspect hydrostatic boundary laws
        c = 0

        % Coefficient multiplying `G` on the \(1/h\) side.
        %
        % - Topic: Inspect hydrostatic boundary laws
        d = 0

        % Coefficient multiplying `G` on the \(h\) side.
        %
        % - Topic: Inspect hydrostatic boundary laws
        e = 0
    end

    methods
        function self = IMHydrostaticBoundaryCondition(options)
            % Create a physical hydrostatic endpoint law.
            %
            % - Topic: Create hydrostatic boundary laws
            % - Declaration: law = IMHydrostaticBoundaryCondition(options)
            % - Parameter options.a: coefficient multiplying `F` on the \(h^0\) side
            % - Parameter options.b: coefficient multiplying `G` on the \(h^0\) side
            % - Parameter options.c: coefficient multiplying `F` on the \(1/h\) side
            % - Parameter options.d: coefficient multiplying `G` on the \(1/h\) side
            % - Parameter options.e: coefficient multiplying `G` on the \(h\) side
            % - Returns law: hydrostatic endpoint law
            arguments
                options.a (1,1) double {mustBeReal, mustBeFinite} = 0
                options.b (1,1) double {mustBeReal, mustBeFinite} = 0
                options.c (1,1) double {mustBeReal, mustBeFinite} = 0
                options.d (1,1) double {mustBeReal, mustBeFinite} = 0
                options.e (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            coefficients = [options.a options.b options.c options.d options.e];
            if all(coefficients == 0)
                error("IMHydrostaticBoundaryCondition:DegenerateLaw", "At least one hydrostatic boundary-law coefficient must be nonzero.");
            end

            self.a = options.a;
            self.b = options.b;
            self.c = options.c;
            self.d = options.d;
            self.e = options.e;
        end

        function boundary = canonicalBoundary(self, options)
            % Convert to a canonical scalar boundary condition.
            %
            % `canonicalBoundary` converts the physical hydrostatic law to
            % the `IMBoundaryCondition` object required by a canonical
            % scalar EVP. Use `formulation="F"` for hydrostatic `F` EVPs
            % and `formulation="G"` for hydrostatic `G` EVPs. The
            % conversion follows the constructor-shaped rules in the class
            % overview.
            %
            % - Topic: Convert hydrostatic boundary laws
            % - Declaration: boundary = canonicalBoundary(law,formulation="G",g=9.81)
            % - Parameter options.formulation: target solved variable, `"F"` or `"G"`
            % - Parameter options.g: gravitational acceleration
            % - Returns boundary: canonical boundary condition
            arguments
                self IMHydrostaticBoundaryCondition
                options.formulation {mustBeTextScalar, mustBeMember(options.formulation, ["F", "G"])}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
            end

            formulation = string(options.formulation);
            g = options.g;
            switch formulation
                case "F"
                    if self.e ~= 0
                        error("IMHydrostaticBoundaryCondition:NonlinearBoundaryLaw", "This hydrostatic boundary law has e=%g, so formulation=""F"" would be nonlinear in lambda=1/h. Set e=0 or convert with formulation=""G"" when d=0.", self.e);
                    end
                    boundary = IMBoundaryCondition(a=-self.a, b=self.b, c=-self.c/g, d=self.d/g);
                case "G"
                    if self.d ~= 0
                        error("IMHydrostaticBoundaryCondition:NonlinearBoundaryLaw", "This hydrostatic boundary law has d=%g, so formulation=""G"" would be nonlinear in lambda=1/h. Set d=0 or convert with formulation=""F"" when e=0.", self.d);
                    end
                    boundary = IMBoundaryCondition(a=self.e, b=self.a, c=self.b/g, d=self.c/g);
            end
        end
    end
end
