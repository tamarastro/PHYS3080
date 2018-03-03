function HzI = Hz_inverse(z,OM,OL) % this is actually Hz/H0
Hz = sqrt((1+z)^2*(OM*z+1)-OL*z*(z+2));
HzI=1.0/Hz;