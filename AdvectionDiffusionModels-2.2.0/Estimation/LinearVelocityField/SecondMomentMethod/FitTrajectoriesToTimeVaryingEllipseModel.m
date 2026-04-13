function [parameters,error,B] = FitTrajectoriesToTimeVaryingEllipseModel( x, y, t, model, dof, S)

[~, ~, q, r] = CenterOfMass( x, y );
[Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);

[parameters,error, B] = FitSecondMomentToTimeVaryingEllipseModel( Mxx, Myy, Mxy, t, model, dof, S);

end
