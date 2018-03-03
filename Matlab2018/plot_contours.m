% plot_contours.m
%
% Description: 
%   Takes the output of fit_models.m
%   And uses it to plot likelihood contours  in the OM-OL plane
%
% Inputs:
%     fit_models_output.mat
%
close all

% Load the data
load fit_models_output.mat

% Calculate the best fit model and print to screen
chi2best = min(min(chi2))
[iOM,iOL]= find(chi2 == chi2best)
OMbest   = OM(iOM)
OLbest   = OL(iOL)

% Define array giving the values of chi2 for 1, 2, and 3 sigma 
% confidence intervals
sig321_2dof = [11.83 6.18 2.296 0.]; % two degrees of freedom
sig321_1dof = [9.    4.   1.    0.]; % one degree of freedom

% Plot contours
clf
hold on
set(gca,'box','on'); % THIS PUTS A BOX AROUND THE PLOT
set(0, 'DefaultAxesFontSize', 14);

[C,h]=contour(OM,OL,chi2'-chi2best,sig321_2dof,'Linewidth',1.2);
colormap winter
%h=axes('FontSize',20)
%set(h,'Fontsize',20)
xlabel('\Omega_M','FontSize',20)
ylabel('\Omega_\Lambda','FontSize',20)
			  
plot(OMbest,OLbest,'x')




