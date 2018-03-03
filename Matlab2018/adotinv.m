function adi = adotinv(a,OM,OL)
%ADOTINV Calculates 1/dot(a)H0
%   The factor of H0 has been divided out.
%   Integrate this function over a to get t in units of Hubble time
%   Multiply by H0 to get the units you prefer (years, seconds...)

OK = 1-OM-OL;
ad = sqrt(OK + OM./a + OL.*a.^2);
adi = 1./ad;

end
