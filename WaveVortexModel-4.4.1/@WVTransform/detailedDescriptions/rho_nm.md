Diagnosed no-motion density profile, `[Nz 1]`, in $$\mathrm{kg\,m^{-3}}$$.

`rho_nm` approximates the statically rearranged current density distribution by fitting its volume-weighted density moments. It is the default reference for `eta_true`, `ape`, and `apv`, with material height $$s$$ defined by $$\rho_{\mathrm{nm}}(s) = \rho_{\mathrm{total}}$$ and displacement $$\eta_{\mathrm{true}} = z - s$$.

A horizontally uniform, strictly stable density field is already its own no-motion profile and is returned exactly. Other states use a bounded damped least-squares fit with strictly ordered interior nodes and the current density extrema as endpoints. The default solver requires no Optimization Toolbox. It rejects failed fits, nonfinite or nonmonotone profiles, and maximum normalized moment residuals above `1e-8` before caching a result. A coarse vertical grid may not resolve a given density distribution to this tolerance; a small optimizer step alone does not certify the result.

A monotone cubic representation supplies both the inverse used by displacement and the integral used by APE. Density plateaus have no unique material height and are rejected. Densities outside the selected profile range are rejected apart from an eight-ulp allowance for density roundoff. APV differentiates the resulting material height through `eta_true`.
