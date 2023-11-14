function [fitresult, gof] = createPoly2Fit(xgrid, x1)
%CREATEFIT(XGRID,X1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: xgrid
%      Y Output: x1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 10-Nov-2023 16:44:04


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xgrid, x1 );

% Set up fittype and options.
ft = fittype( 'poly2' );
excludedPoints = (xData < 0.007) | (xData > 0.01);
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'x1 vs. xgrid', 'Excluded x1 vs. xgrid', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'xgrid', 'Interpreter', 'none' );
ylabel( 'x1', 'Interpreter', 'none' );
grid on


