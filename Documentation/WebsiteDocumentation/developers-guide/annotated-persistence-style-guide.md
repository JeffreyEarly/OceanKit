---
layout: default
title: Annotated persistence
parent: Developers guide
nav_order: 4
---

# Annotated persistence style guide

This guide defines the shared OceanKit convention for classes that persist state through `CAAnnotatedClass` and NetCDF-backed files.

The goal is not to force every package to use identical scientific types or identical constructor signatures. The goal is to make annotated persistence predictable enough that package authors, tests, documentation tooling, and AI-assisted edits can move between repositories without relearning the structure each time.

The shared convention should be as simple as the class allows. `wave-vortex-model` uses the most elaborate form of this pattern because many of its constructors reconstruct positional scientific identities from several composed roles. Most packages should use a lighter form.

For simple persisted classes, OceanKit prefers `name=value` constructors whose argument names match `classRequiredPropertyNames()`. That alignment keeps the persistence path small, direct, and easier to maintain.

When a class has both a cheap canonical persisted state and an expensive scientific setup path, OceanKit prefers making the constructor the cheap canonical constructor and moving the expensive setup into explicit static factories that delegate to that constructor.

Use this guide together with the MATLAB style guide and documentation style guide. Package-specific documentation may extend this pattern, but it should not redefine the generic contract differently.

## Choose the lightest pattern that fits

Use the simplest persistence pattern that preserves a clear constructor contract:

1. Direct annotated reconstruction.
   Use this when the constructor can be called directly from the persisted required-property set. This is the preferred pattern for simple classes, and it works best when the constructor uses `name=value` arguments whose names match `classRequiredPropertyNames()`. In this case `CAAnnotatedClass.annotatedClassFromFile()` or `annotatedClassFromGroup()` is sufficient, and a thin class-named wrapper is optional.

2. Thin adapter reconstruction.
   Use this when the persisted properties are simple, but the constructor expects positional inputs, canonicalized options, or a small amount of argument shaping. In this case keep `classDefinedPropertyAnnotations()` and `classRequiredPropertyNames()` on the class, and add a small explicit `fromGroup(...)` adapter that turns persisted values into constructor inputs.

3. Role-based reconstruction.
   Use this when persistence is composed across several scientific roles, parents, or subsystems. This is the `wave-vortex-model` pattern with helpers such as `propertyAnnotationsFor<Role>()`, `namesOfRequiredPropertiesFor<Role>()`, and `requiredPropertiesFor<Role>FromGroup(...)`.

Default to level 1 or 2 unless the class genuinely needs level 3.

## Required contract

Use annotated persistence when a class writes its state to NetCDF through `CAAnnotatedClass` metadata and is expected to be reconstructed later.

For each public persisted object type:

- subclass `CAAnnotatedClass`
- implement `classDefinedPropertyAnnotations()`
- implement `classRequiredPropertyNames()`
- provide a public instance `writeToFile(path, ...)`, either inherited or overridden

Add explicit `roleFromFile(path, options)` and `roleFromGroup(group)` factories when one of the following is true:

- the constructor cannot be called directly from the required-property set
- the class needs package-specific validation or post-construction restoration
- the package wants a stable role-specific public API such as `geometryFromFile(...)` or `distributionFromGroup(...)`

The instance-side writer and the class-side reader are a pair. `writeToFile(...)` saves the current object state. `roleFromFile(...)` and `roleFromGroup(...)` reconstruct a new object from persisted state.

Avoid hidden file I/O in constructors. Constructors should keep their normal scientific initialization contract. File reconstruction belongs in explicit factories.

## Default pattern for simple classes

For simple persisted classes, keep the structure local to the class:

- a constructor that prefers `name=value` arguments matching `classRequiredPropertyNames()`
- when the class also has an expensive scientific setup path, explicit static factories should build that state and then delegate to the constructor
- `classDefinedPropertyAnnotations()`
- `classRequiredPropertyNames()`
- optional thin `classNameFromFile(path)` wrapper
- optional thin `classNameFromGroup(group)` adapter

In this pattern, the class annotations and required-property list are the primary persistence contract. Do not introduce extra role helper methods unless they actually clarify reconstruction.

Treat this constructor alignment as the default OceanKit design for simple persisted classes. It keeps constructor calls, persistence metadata, generic annotated reconstruction, tests, and documentation all speaking the same vocabulary.

For classes whose ordinary scientific setup is materially more expensive than restart, prefer a public canonical constructor that takes the persisted-state vocabulary directly, and explicit static factories such as `fromGriddedValues(...)` or `fromKnotPoints(...)` for the scientific setup. Let those factories validate their own inputs, do the expensive work, and then delegate to the constructor.

Keep formal `arguments` validation in place. Cheap constructors should validate canonical state directly, and scientific factories should validate their own source-specific inputs rather than weakening either contract.

If a class has a strong scientific reason to keep positional constructor inputs, that is still allowed. In that case, use the thin-adapter pattern explicitly rather than forcing the persistence layer to guess how stored fields map back onto positional arguments.

Avoid package-private orchestration helpers with names like `reconstructPersistedRole(...)` or `readConcreteRoleFromFile(...)` unless the logic is genuinely shared and would otherwise be duplicated in a way that hurts readability. In most packages, the explicit `roleFromFile(...)` and `roleFromGroup(...)` methods are clearer.

Prefer short inline constructor branching over parser helpers such as `parseBootstrapInputs(...)` or `shouldUseBootstrapConstruction(...)` when the logic is truly small. Once a class starts carrying distinct scientific-setup and canonical-restart paths, prefer explicit factories over growing constructor mode logic.

## Helper naming for advanced cases

When a class needs role-based reconstruction, use the following helper names for each persisted role:

- `propertyAnnotationsFor<Role>()`
- `namesOfRequiredPropertiesFor<Role>()`
- `requiredPropertiesFor<Role>FromGroup(group, options)`

These helpers separate three concerns:

- `propertyAnnotationsFor<Role>()` declares the persisted dimensions, variables, objects, functions, and attributes
- `namesOfRequiredPropertiesFor<Role>()` names the saved properties required to recover equivalent object state
- `requiredPropertiesFor<Role>FromGroup(...)` adapts stored values into the constructor shape used by the class, such as `[Lxy, Nxy, options]`

`classDefinedPropertyAnnotations()` and `classRequiredPropertyNames()` should usually be thin wrappers around the role-specific helpers in this advanced pattern.

For composed subclasses, the following optional helper names are recommended when they improve clarity:

- `newRequiredPropertyNames()`
- `newNonrequiredPropertyNames()`

Use these to make subclass additions and intentional removals explicit rather than burying them in a long `union(...)` or `setdiff(...)` expression.

## Required properties and reconstruction

`classRequiredPropertyNames()` is a persistence contract, not a copy of the constructor signature.

List the persisted properties needed to recover equivalent object state, including values that may be optional during ordinary construction but are required to reconstruct the saved state exactly.

For simple classes, the preferred design is for these required-property names to also be the constructor's `name=value` argument names. When that is true, the generic annotated read path usually stays small enough that no additional reconstruction adapter is needed.

If only a tiny amount of extra bootstrap state is needed, additional validated `options.<name>` fields are acceptable. Once that extra state creates a second substantive initialization path, prefer the cheap canonical constructor plus explicit static factories.

In the direct and thin-adapter patterns, `roleFromGroup(...)` itself may be the adapter between persisted state and constructor inputs.

In the role-based pattern, `requiredPropertiesFor<Role>FromGroup(...)` is the adapter between persisted state and constructor inputs. This helper may:

- validate that the required persisted properties are present
- read the stored values from a `NetCDFGroup`
- derive constructor sizes or positional inputs from saved arrays
- split the result into positional constructor inputs plus name-value options

This split is important when constructors use a compact public API such as required positional arguments followed by name-value options, but the file format stores a larger set of explicit state variables.

## Read path

Use the following read flow when explicit factories are needed:

1. `roleFromFile(path, options)` opens the NetCDF file and delegates to `roleFromGroup(group)`.
2. `roleFromGroup(group)` validates the persisted state with `CAAnnotatedClass.throwErrorIfMissingProperties(...)` or `CAAnnotatedClass.canInitializeDirectlyFromGroup(...)`.
3. `roleFromGroup(group)` calls `requiredPropertiesFor<Role>FromGroup(...)` to fetch and adapt the persisted values.
4. `roleFromGroup(group)` constructs the object with the normal class constructor.

Use `CAAnnotatedClass.propertyValuesFromGroup(...)` to read persisted values rather than reimplementing NetCDF lookup logic in each package.

The generic factories `CAAnnotatedClass.annotatedClassFromFile()` and `CAAnnotatedClass.annotatedClassFromGroup()` are appropriate when the constructor can be called directly from the required-property set without package-specific assembly. If reconstruction needs argument shaping, derived positional inputs, nested polymorphic dispatch, or additional initialization, expose explicit `roleFromFile(...)` and `roleFromGroup(...)` factories instead.

## Write path

`writeToFile(path, ...)` is the public instance-side persistence entry point.

By default, the writer should rely on the annotated-property contract:

- determine the properties to persist
- map them to `CAPropertyAnnotation` instances
- delegate the actual NetCDF write through `CAAnnotatedClass`

Subclasses may override `writeToFile(...)` or `writeToGroup(...)` when they need to add domain-specific variables, groups, or metadata. Those overrides should still build on the shared annotated persistence path rather than inventing a separate serialization mechanism for the same object.

Keep file options explicit and stable. Prefer names such as `shouldOverwriteExisting`, `shouldReadOnly`, and `iTime`.

## Inheritance and composition

Root classes should own their role helpers when the advanced role-based pattern is in use. In particular, the root class for a persisted role should define:

- `propertyAnnotationsFor<Role>()`
- `namesOfRequiredPropertiesFor<Role>()`
- `requiredPropertiesFor<Role>FromGroup(...)`
- `roleFromFile(...)`
- `roleFromGroup(...)`

Subclasses should combine parent persistence contracts explicitly:

- combine required-property sets with `union(...)`
- remove intentionally forced or derived values with `setdiff(...)`
- combine annotation arrays explicitly with `cat(...)`

Do not rely on implicit inheritance of persistence structure. Make the combined contract visible in code so a reader can tell which parent contributes which persisted state.

When a subclass makes an assumption that replaces a parent property, remove that property from the required set rather than persisting it as dead state.

## Templates

### Simple persisted class

```matlab
classdef ExampleDistribution < CAAnnotatedClass
    properties (SetAccess = private)
        sigma
        nu
    end

    methods
        function self = ExampleDistribution(options)
            arguments
                options.sigma (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 1
                options.nu (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive} = 4
            end
            self.sigma = options.sigma;
            self.nu = options.nu;
        end
    end

    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('sigma', {}, '', 'scale parameter');
            propertyAnnotations(end+1) = CANumericProperty('nu', {}, '', 'degrees of freedom');
        end

        function names = classRequiredPropertyNames()
            names = {'sigma', 'nu'};
        end
    end
end
```

### Cheap canonical constructor plus scientific factory

```matlab
classdef ExampleSpline < CAAnnotatedClass
    properties (SetAccess = private)
        gridAxes
        knotAxes
        xi
    end

    methods
        function self = ExampleSpline(options)
            arguments
                options.gridAxes (:,1) ExampleAxis
                options.knotAxes (:,1) ExampleAxis
                options.xi {mustBeNumeric,mustBeReal,mustBeFinite}
            end
            self.gridAxes = options.gridAxes;
            self.knotAxes = options.knotAxes;
            self.xi = options.xi;
        end
    end

    methods (Static)
        function self = fromGriddedValues(gridVectors, values, options)
            arguments
                gridVectors
                values {mustBeNumeric,mustBeReal,mustBeFinite}
                options.S {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative} = 3
            end

            [gridAxes, knotAxes, xi] = ExampleSpline.solveFromValues(gridVectors, values, options.S);
            self = ExampleSpline(gridAxes=gridAxes, knotAxes=knotAxes, xi=xi);
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAObjectProperty('gridAxes', 'Grid axes.');
            propertyAnnotations(end+1) = CAObjectProperty('knotAxes', 'Knot axes.');
            propertyAnnotations(end+1) = CANumericProperty('xi', {}, '', 'Coefficient state.');
        end

        function names = classRequiredPropertyNames()
            names = {'gridAxes', 'knotAxes', 'xi'};
        end
    end
end
```

### Composed or multiple-inheritance class

```matlab
classdef ExampleBarotropicGeometry < ExampleGeometry
    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = ExampleBarotropicGeometry.propertyAnnotationsForGeometry();
        end

        function names = classRequiredPropertyNames()
            names = ExampleBarotropicGeometry.namesOfRequiredPropertiesForGeometry();
        end

        function names = namesOfRequiredPropertiesForGeometry()
            names = ExampleGeometry.namesOfRequiredPropertiesForGeometry();
            names = union(names, ExampleRotation.namesOfRequiredPropertiesForRotation());
            names = union(names, ExampleBarotropicGeometry.newRequiredPropertyNames());
            names = setdiff(names, ExampleBarotropicGeometry.newNonrequiredPropertyNames());
        end

        function names = newRequiredPropertyNames()
            names = {'h', 'j'};
        end

        function names = newNonrequiredPropertyNames()
            names = {'Nz'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = ExampleGeometry.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2, propertyAnnotations, ExampleRotation.propertyAnnotationsForRotation());
            propertyAnnotations(end+1) = CANumericProperty('h', {}, 'm', 'equivalent depth');
            propertyAnnotations(end+1) = CANumericProperty('j', {}, '', 'mode number');
        end

        function [Lxy, Nxy, options] = requiredPropertiesForGeometryFromGroup(group)
            [Lxy, Nxy, geomOptions] = ExampleGeometry.requiredPropertiesForGeometryFromGroup(group);
            rotationOptions = ExampleRotation.requiredPropertiesForRotationFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, ExampleBarotropicGeometry.newRequiredPropertyNames());
            newOptions = namedargs2cell(vars);
            options = cat(2, geomOptions, rotationOptions, newOptions);
        end
    end
end
```

## Existing reference patterns

Reference implementations for the three levels live in:

- `class-annotations/@CAAnnotatedClass`
- simple or thin-adapter package-local classes such as `distributions/@NormalDistribution` and `distributions/@StudentTDistribution`
- `wave-vortex-model/@WVGeometryDoublyPeriodic`
- `wave-vortex-model/@WVGeometryDoublyPeriodicBarotropic`
- `wave-vortex-model/@WVTransformBarotropicQG`

Use those classes as concrete examples when the abstract rule needs a repository-backed reference.
