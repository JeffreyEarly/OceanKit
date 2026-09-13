Set a barotropic geostrophic state from sea-surface height.

`ssh` is a function handle evaluated on the horizontal grid. The transform converts height to geostrophic streamfunction using $$\psi=g\eta/f$$. Set `shouldRemoveMeanPressure=true` to remove the horizontal mean first; the default is `false`.

```matlab
wvt.setSSH(@(x,y)0.1*cos(2*pi*x/wvt.Lx));
```

- Declaration: setSSH(ssh,options)
- Parameter ssh: function handle accepting horizontal coordinate arrays and returning height in meters
- Parameter options.shouldRemoveMeanPressure: remove the horizontal mean before initialization; default `false`
