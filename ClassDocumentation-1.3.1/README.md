# ClassDocumentation

ClassDocumentation generates a Markdown class reference from MATLAB class metadata and structured help comments. It writes a class landing page and one page for each documented public method or property.

```matlab
documentation = ClassDocumentation("ExampleClass", ...
    buildFolder="docs", ...
    websiteFolder="classes", ...
    parent="Transforms", ...
    grandparent="Class documentation", ...
    methodGrandparent="Transforms");
documentation.writeToFile();
```

`parent` and `grandparent` configure the class landing page. Each generated method or property page uses the class name as its immediate parent and `methodGrandparent` as its grandparent. `methodGrandparent` defaults to `"Classes"` so existing documentation builds retain their current navigation until they opt into a different hierarchy.

For annotated properties, ClassDocumentation preserves the original machine-readable `units` value and formats recognized CF/UDUNITS-style expressions for display. For example, `kg m-1 s-2` renders as upright, spaced MathJax; legacy spellings such as `kg/m/s2` remain supported during migration. Dimensionless and directional units receive readable prose, while unrecognized values are shown literally rather than reinterpreted.
