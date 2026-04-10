---
layout: default
title: Annotated persistence
parent: Developers guide
nav_order: 4
---

# Annotated persistence style guide

This guide defines the shared OceanKit convention for classes that persist state through `CAAnnotatedClass` and NetCDF-backed files.

The goal is not to force every package to use identical scientific types or identical constructor signatures. The goal is to make annotated persistence predictable enough that package authors, tests, documentation tooling, and AI-assisted edits can move between repositories without relearning the structure each time.

Use this guide together with the MATLAB style guide and documentation style guide. Package-specific documentation may extend this pattern, but it should not redefine the generic contract differently.

## Required contract

Use annotated persistence when a class writes its state to NetCDF through `CAAnnotatedClass` metadata and is expected to be reconstructed later.

For each public persisted object type:

- subclass `CAAnnotatedClass`
- implement `classDefinedPropertyAnnotations()`
- implement `classRequiredPropertyNames()`
- provide an instance `writeToFile(path, ...)`
- provide an explicit source-specific `roleFromFile(path, options)` factory
- provide an explicit `roleFromGroup(group)` factory

The instance-side writer and the class-side reader are a pair. `writeToFile(...)` saves the current object state. `roleFromFile(...)` and `roleFromGroup(...)` reconstruct a new object from persisted state.

Avoid hidden file I/O in constructors. Constructors should keep their normal scientific initialization contract. File reconstruction belongs in explicit factories.

## Helper naming

Use the following helper names for each persisted role:

- `propertyAnnotationsFor<Role>()`
- `namesOfRequiredPropertiesFor<Role>()`
- `requiredPropertiesFor<Role>FromGroup(group, options)`

These helpers separate three concerns:

- `propertyAnnotationsFor<Role>()` declares the persisted dimensions, variables, objects, functions, and attributes
- `namesOfRequiredPropertiesFor<Role>()` names the saved properties required to recover equivalent object state
- `requiredPropertiesFor<Role>FromGroup(...)` adapts stored values into the constructor shape used by the class, such as `[Lxy, Nxy, options]`

`classDefinedPropertyAnnotations()` and `classRequiredPropertyNames()` should usually be thin wrappers around the role-specific helpers.

For composed subclasses, the following optional helper names are recommended when they improve clarity:

- `newRequiredPropertyNames()`
- `newNonrequiredPropertyNames()`

Use these to make subclass additions and intentional removals explicit rather than burying them in a long `union(...)` or `setdiff(...)` expression.

## Required properties and reconstruction

`classRequiredPropertyNames()` is a persistence contract, not a copy of the constructor signature.

List the persisted properties needed to recover equivalent object state, including values that may be optional during ordinary construction but are required to reconstruct the saved state exactly.

`requiredPropertiesFor<Role>FromGroup(...)` is the adapter between persisted state and constructor inputs. This helper may:

- validate that the required persisted properties are present
- read the stored values from a `NetCDFGroup`
- derive constructor sizes or positional inputs from saved arrays
- split the result into positional constructor inputs plus name-value options

This split is important when constructors use a compact public API such as required positional arguments followed by name-value options, but the file format stores a larger set of explicit state variables.

## Read path

Use the following read flow for explicit factories:

1. `roleFromFile(path, options)` opens the NetCDF file and delegates to `roleFromGroup(group)`.
2. `roleFromGroup(group)` validates the persisted state with `CAAnnotatedClass.throwErrorIfMissingProperties(...)` or `CAAnnotatedClass.canInitializeDirectlyFromGroup(...)`.
3. `roleFromGroup(group)` calls `requiredPropertiesFor<Role>FromGroup(...)` to fetch and adapt the persisted values.
4. `roleFromGroup(group)` constructs the object with the normal class constructor.

Use `CAAnnotatedClass.propertyValuesFromGroup(...)` to read persisted values rather than reimplementing NetCDF lookup logic in each package.

The generic factories `CAAnnotatedClass.annotatedClassFromFile()` and `CAAnnotatedClass.annotatedClassFromGroup()` are only appropriate when the constructor can be called directly from the required-property set without package-specific assembly. If reconstruction needs role-specific argument shaping, derived positional inputs, or additional initialization, expose explicit `roleFromFile(...)` and `roleFromGroup(...)` factories instead.

## Write path

`writeToFile(path, ...)` is the public instance-side persistence entry point.

By default, the writer should rely on the annotated-property contract:

- determine the properties to persist
- map them to `CAPropertyAnnotation` instances
- delegate the actual NetCDF write through `CAAnnotatedClass`

Subclasses may override `writeToFile(...)` or `writeToGroup(...)` when they need to add domain-specific variables, groups, or metadata. Those overrides should still build on the shared annotated persistence path rather than inventing a separate serialization mechanism for the same object.

Keep file options explicit and stable. Prefer names such as `shouldOverwriteExisting`, `shouldReadOnly`, and `iTime`.

## Inheritance and composition

Root classes should own their role helpers. In particular, the root class for a persisted role should define:

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
classdef ExampleGeometry < CAAnnotatedClass
    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = ExampleGeometry.propertyAnnotationsForGeometry();
        end

        function names = classRequiredPropertyNames()
            names = ExampleGeometry.namesOfRequiredPropertiesForGeometry();
        end

        function names = namesOfRequiredPropertiesForGeometry()
            names = {'x', 'y', 'Lx', 'Ly', 'shouldAntialias'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('x', 'm', 'x coordinate');
            propertyAnnotations(end+1) = CADimensionProperty('y', 'm', 'y coordinate');
            propertyAnnotations(end+1) = CANumericProperty('Lx', {}, 'm', 'domain length in x');
            propertyAnnotations(end+1) = CANumericProperty('Ly', {}, 'm', 'domain length in y');
        end

        function [Lxy, Nxy, options] = requiredPropertiesForGeometryFromGroup(group)
            vars = CAAnnotatedClass.propertyValuesFromGroup( ...
                group, ExampleGeometry.namesOfRequiredPropertiesForGeometry());
            Nxy = [length(vars.x), length(vars.y)];
            Lxy = [vars.Lx, vars.Ly];
            vars = rmfield(vars, {'x', 'y', 'Lx', 'Ly'});
            options = namedargs2cell(vars);
        end

        function geometry = geometryFromFile(path)
            ncfile = NetCDFFile(path);
            geometry = ExampleGeometry.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            CAAnnotatedClass.throwErrorIfMissingProperties( ...
                group, ExampleGeometry.namesOfRequiredPropertiesForGeometry());
            [Lxy, Nxy, options] = ExampleGeometry.requiredPropertiesForGeometryFromGroup(group);
            geometry = ExampleGeometry(Lxy, Nxy, options{:});
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
            vars = CAAnnotatedClass.propertyValuesFromGroup( ...
                group, ExampleBarotropicGeometry.newRequiredPropertyNames());
            newOptions = namedargs2cell(vars);
            options = cat(2, geomOptions, rotationOptions, newOptions);
        end
    end
end
```

## Existing reference patterns

The current reference implementations for this style live in:

- `class-annotations/@CAAnnotatedClass`
- `wave-vortex-model/@WVGeometryDoublyPeriodic`
- `wave-vortex-model/@WVGeometryDoublyPeriodicBarotropic`
- `wave-vortex-model/@WVTransformBarotropicQG`

Use those classes as concrete examples when the abstract rule needs a repository-backed reference.
