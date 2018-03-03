function DM = dist_mod(zz,OM,OL)
% This actually returns the Hubble-constant free distance modulus.
% i.e. the actual distance modulus is DM + 5*Alog10(c/H0) 
%
%
% Note that the input z can be an array.  
% The output is then an array of distance moduli.
% Functions called:
%   Hz_inverse(z,OM,OL):
%       this is just 1/H(z)

% Define the necessary constants and functions.
% OK is OmegaK, the curvature
% x is comoving distance
% D is proper distance

    OK = 1.0 - OM - OL;
    
    % Integrate H(z) over z to get comoving distance, x
    for ii=1:numel(zz)
        z     = zz(ii);
        x(ii) = quadv(@(z)Hz_inverse(z,OM,OL),0,z);
    end
    x=x';  % ' gives transpose of vector
    
    % Convert this to proper distance 
    if(OK < 0.0)
        R0 = 1/sqrt(-OK);
        D = R0.*sin(x./R0);
    end
    if(OK > 0.0)
        R0 = 1/sqrt(OK);
        D = R0.*sinh(x./R0);
    end
    if(OK == 0.0)
        D=x;
    end

    % Convert proper distance to luminosity distance
    LumDist = D.*(1+zz);
    
    % Calculate distance modulus from luminosity distance
    DM = 5.*log10(LumDist);