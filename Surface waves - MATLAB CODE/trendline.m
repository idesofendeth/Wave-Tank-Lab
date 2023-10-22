function [px,py] = trendline(x,err)

p = polyfit(x, err, 1);
px = [min(x) max(x)];
py = polyval(p, px);

end