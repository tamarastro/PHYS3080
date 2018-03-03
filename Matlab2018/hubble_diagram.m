% hubble_diagram.m
%
% Description: 
%   Plot a magnitude-redshift diagram.
%
% Inputs:
%   jla_lcparams_simple_new.txt  (Results from Joint Light Curve Analysis, Betoule et al. 2014)
%
% Outputs:
%
% User defined functions called:
%   dist_mod(z,OM,OL): 
%       calculates distance modulus with parameters z, OM, and OL
%       needs Hz_inverse.m
% 
% MATLAB functions called:
%


%-------------------
% Read in data
%-------------------
clear
fid=fopen('jla_lcparams_simple_new.txt');
data=textscan(fid,'%f%f%f%f','headerlines',15);
fclose(fid);

% Supernova details
zz    =data{1};  % Redshift
mu    =data{2};  % Distance modulus - abs mag
muerr =data{3};  % Distance modulus uncertainty
id    =data{4}; % Source ID 
                 %  [3 = Nearby sample           ]
                 %  [2 = Sloan Digital Sky Survey]
                 %  [1 = Supernova Legacy Survey ]
                 %  [4 = Hubble Space Telescope  ]
id_nums  = [3 2 1 4];
id_names = ['Nearby','SDSS','SNLS','HST'];

%---------------------------------
% Calculate a model prediction
%---------------------------------
% Define an array of redshifts to calculate the model for
zz_model = logspace(-2,0.2,50)'; %Make it logarithmic to better sample low-z
% Calculate the distance modulus for those redshifts  
mu_model=dist_mod(zz_model,0.3,0.7);

% To calculate the arbitrary offset we could use:
% total((mu_data-mu_theory)/muerr^2)/total(1./muerr^2)
% But for now I'm just going to set it manually
mscr=43.1 ;
mu_model = mu_model+mscr ;

%-------------------
% PLOT DATA WITH ERROR BARS
%-------------------
% Clear existing graphs and create axes
clf
axes1 = axes('FontSize',16);
box('on');
hold('all'); 

% Plot the data (different data sources in different colours)
for i=[1:4] 
    index = find(id == id_nums(i));
    errorbar(zz(index),mu(index),muerr(index),'.')
end
legend('Nearby','SDSS','SNLS','HST','location','SouthEast')
legend('boxoff')

plot(zz_model,mu_model,'-k','Linewidth',1.2)
% Add annotations to plot 
xlabel('Redshift' ,'FontSize',16);
ylabel('Peak magnitude','FontSize',16);
xlim  ([0.0 1.4])
% Legend
text(0.1,36,'Model: (\Omega_M,\Omega_\Lambda)=(0.3,0.7)','FontSize',16,'FontName','Times New Roman')
text(0.1,35,'Data: SDSS compilation 2009','FontSize',16,'FontName','Times New Roman')

hold('off');
