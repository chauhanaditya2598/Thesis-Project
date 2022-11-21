%% Let's create function to calculate the buckling ratio (B)
% Ef = Young's modulus of fiber
% df = diameter of fiber
% zeta = hydrodynamic drag coefficient
% gamma = shear rate
% eta_m = viscosity
% l = current fiber length

function [B,Lub] = BuckRatio(Ef,df,zeta,eta_m,gamma,l)

Lub = (((pi^3).*Ef.*(df)^4)/(4*zeta.*eta_m.*gamma))^(1/4);

B = ((l./Lub).^4);

end
