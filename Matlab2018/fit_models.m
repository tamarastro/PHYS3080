% fit_models.m
%
% Description: 
%   Find the best-fit model to the supernova data 
%   The supernova data are distance moduli and redshifts from SDSS09
%
% Inputs:
%   jla_lcparams_simple.txt  (Results from Joint Light Curve Analysis, Betoule et al. 2014)
%
% Outputs:
%
% User defined functions called:
%   dist_mod(z,OM,OL): 
%       calculates distance modulus for z, OM, and OL
%       z can be an array
%
%   Hz_inverse(z,OM,OL):
%       1/H(z)... this gets called by dist_mod
% 
%
'Running fit_models.m'

%-------------------
% Read in data
%-------------------
clear
fid=fopen('jla_lcparams_simple_new.txt');
data=textscan(fid,'%f%f%f%f','headerlines',15);
fclose(fid);

% Supernova details
zz    =data{1};  % Redshift
zzn   =numel(zz) % Count number of redshifts
mu    =data{2};  % Distance modulus - abs mag
muerr =data{3};  % Distance modulus uncertainty
muerr = muerr %+ 0.1; % Add intrinsic diversity in SNe that is not taken into account in the uncertainties provided
id    =data{4}; % Source ID 
                 %  [3 = Nearby sample           ]
                 %  [2 = Sloan Digital Sky Survey]
                 %  [1 = Supernova Legacy Survey ]
                 %  [4 = Hubble Space Telescope  ]
id_nums  = [3 2 1 4];
id_names = ['Nearby','SDSS','SNLS','HST'];


%--------------------------------------------------
% Define vectors for the parameters we are testing
%--------------------------------------------------
% Note that the integration in this simple program fails for bounce models 
% (models with no beginning in time) so either avoid OLmax > 1.15 or so, 
% or re-write the program to make it more robust. 
% The numbers are fairly small at the moment.  Increase the resolution to 
% get better fits and prettier contours.
% Initialise values of OM and OL to loop over
OM      = [0.0 : 0.05 : 0.7]';
OMn     = numel(OM);
OL      = [0.0 : 0.05 : 1.0]';
OLn     = numel(OL);
% Initialise a range of normalisation values to test
mscr_guess = 43.1 %Approx 5.*log10(c_H0)+25-19;
mscr_lo    = mscr_guess-0.5
mscr_hi    = mscr_guess+0.5
mscr       = [mscr_lo : 0.005 : mscr_hi]';
mscrn      = numel(mscr);

% It is going to make the code faster if we calculate mu for a sparse array
% and interpolate, so make an array of redshifts for that calculation.
% Make it logarithmic to better sample low-z for more accurate interpoloation where the curve is curviest!
zz_model = logspace(-2,0.2,50)';
% Initialise a chi2 matrix, and make all values huge so that 
% our first estimate of chi2 is guaranteed to be an improvement.
chi2      = ones (OMn,OLn)*1e20;
% Also introduce a matrix in which to save the normalisation factors.
mscr_used = zeros(OMn,OLn);
% Finally introduce a matrix in which to save the distance moduli.
mu_model  = zeros(zzn);

%-----------------------------------------------------------------------
% Compute probability (chi2) for each combination of OM and OL
%-----------------------------------------------------------------------
'StartingLoops'
for i=[1:OMn] 
    for j=[1:OLn]
%       Calculate distance moduli predicted by the model
%       for the range of z's we need
        mu_model_sparse = dist_mod(zz_model,OM(i),OL(j));
%       Interpolate to the redshifts of the data
        mu_model = interp1(zz_model,mu_model_sparse,zz);
%       FYI, full calculation without interpolation would be much slower:
        %mu_model = dist_mod(zz,OM(i),OL(j));
        
%       Compare this distance modulus to the data to get 
%       the chisquared value.  For each distance modulus
%       we have to check a bunch of normalisations, because
%       we don't know the absolute magnitude of the supenovae
%       well enough, nor the Hubble constant well enough.  
%       We do know the relative magnitudes, which means
%       we fit only for the slope, and not the y-intercept.
%       (Note, this is a very clunky way to test mscr values, 
%       but it makes it obvious what is happening.)
        for k=[1:mscrn]
            % Apply the first normalisation estimate
            mu_model_norm = mu_model + mscr(k);
            % Calculate chi2 for that normalisation
            chi2_test = sum( (mu_model_norm - mu).^2 ./ muerr.^2 );
            %mean(mu_model_norm - mu)
            %pause
            % If that's an improvement, keep it, if not move on.
            if (chi2_test < chi2(i,j)) 
                chi2(i,j)=chi2_test;
                mscr_used(i,j)=mscr(k);
            end
        end
    end
    PercentComplete = (100*i)/OMn
end

% Print out the values for the best fit parameters:
chi2best = min(min(chi2));
bestchi2perdof = chi2best/(zzn-3)
[iOM,iOL]= find(chi2 == chi2best);
OMbest   = OM(iOM)
OLbest   = OL(iOL)
save fit_models_output.mat %zz OM OL chi2 chi2best OMbestfit OLbestfit

% Check that the mscr range we used was adequate
% The values used must not be at the edge of those tested
min_mscr_used=min(min(mscr_used));
max_mscr_used=max(max(mscr_used));
mscr_range_tested = [mscr_lo,mscr_hi];

if (min(min(mscr_used)) <= mscr_lo )
    'WARNING: minimum mscr reached, extend range '
end
if (max(max(mscr_used)) >= mscr_hi )
    'WARNING: maximum mscr reached, extend range'
end

% Repeat the printing, but this time to a file for later reference
fileID=fopen('fit_models_output.txt','w');
fprintf(fileID,'Best Chi2/dof = %6.3f\n',bestchi2perdof);
fprintf(fileID,'Best(OM,OL) = (%6.3f,%6.3f)\n',OMbest,OLbest);
fprintf(fileID,'Mscr range tested: %6.3f < mscr < %6.3f\n',mscr_lo,mscr_hi);
fprintf(fileID,'Mscr range used  : %6.3f < mscr < %6.3f\n',min(min(mscr_used)),max(max(mscr_used)));
fprintf(fileID,'(Mscr used must be within the testing range)');
if (min(min(mscr_used)) <= mscr_lo ) 
    fprintf(fileID,'WARNING: minimum mscr reached (%6.3f), extend range\n',min(min(mscr_used)));
end
if (max(max(mscr_used)) >= mscr_hi )
    fprintf(fileID,'WARNING: maximum mscr reached (%6.3f), extend range\n',max(max(mscr_used)));
end
fclose(fileID);

   
