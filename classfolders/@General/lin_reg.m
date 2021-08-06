function [p_0, p_1, yfit, rsq] = lin_reg(x_vals,y_vals)
p = polyfit(x_vals,y_vals,1);
p_0 = p(1);p_1 = p(2); % y_fit = p_0*x + p_1
yfit = polyval(p,x_vals);

yresid = y_vals - yfit;

SSresid = sum(yresid.^2);

SStotal = (length(y_vals)-1) * var(y_vals);

rsq = 1 - SSresid/SStotal;

end