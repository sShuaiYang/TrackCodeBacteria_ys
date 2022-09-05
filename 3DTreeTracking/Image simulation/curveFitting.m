function b=curveFitting(xdata,ydata)

% Set up fittype and options.
ft = fittype( 'exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = -Inf;
opts.StartPoint = 0.970592781760616;
opts.Upper = Inf;

% Fit model to data.
fitresult= fit(xdata,ydata,ft,opts);
b=fitresult.b;
end
