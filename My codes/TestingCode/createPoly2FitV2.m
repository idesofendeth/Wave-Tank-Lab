function [fitresult, gof] = createPoly2FitV2(xgrid, x1,peakPos,deltax)
%CREATEFIT(XGRID,X1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: xgrid
%      Y Output: x1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%       in fitresult, b is the x position of the fitted peak crest!
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 14-Nov-2023 11:37:05


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xgrid, x1 );

% Set up fittype and options.
ft = fittype( 'a*(x-b)^2+c', 'independent', 'x', 'dependent', 'y' );
excludedPoints = (xData < peakPos-deltax) | (xData > peakPos+deltax);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-9*10^5 0 -Inf];
opts.StartPoint = [0.743132468124916 peakPos 1];
opts.Upper = [0 Inf Inf];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'x1 vs. xgrid', 'Excluded x1 vs. xgrid', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'xgrid', 'Interpreter', 'none' );
% ylabel( 'x1', 'Interpreter', 'none' );
% ylim([0 1])
% xlim([0 xgrid(end)])
% grid on

peakPos
peakPos-deltax
peakPos+deltax


