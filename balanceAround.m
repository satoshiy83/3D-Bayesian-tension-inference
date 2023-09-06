% Function to return coefficients to balance tension and pressures around a
% vertex.
% Written by Satoshi Yamashita.

function [result1,result2] = balanceAround(center,vertices)
% Function to return coefficients to balance tension and pressures.
% [result1,result2] = balanceAround(center,vertices)
% Argument center is a row vector representing a vertex at the center.
% Argument vertices is an array of vertices adjacent to the center.
% Return value result1 is a matrix whose first row represent coefficients
% for the tensions in one arbitrary direction, and second row represent
% coefficients for the tensions in the orthogonal direction.
% Return value result2 is a matrix whose first row represent coefficients
% for the pressures in the first direction, and second row represent
% coefficients for the pressures in the second direction.
% Let edge between center and vertices(i) be denoted by e_i and tension
% exerted on it be denoted by T_i. Then an i-th column of result1
% represents the coefficient for T_i, and i-th column of result2 represents
% the coefficient for a cell which contains e_i and e_{i + 1}.

angles = angleAround(center,vertices);
lengths = edgeLengthAround(center,vertices);

%% Incident angles
mengths = sqrt(sum((center - vertices) .^ 2,2));
correct = lengths ./ mengths;
straight_edge = true;

%% Coefficients for tensions
bngles = cumsum(angles);
bngles = [0; bngles(1:end - 1)];

coef_x = -cos(bngles);
coef_y = -sin(bngles);

if straight_edge
    result1 = [coef_x'; coef_y'] .* correct';
else
    result1 = [coef_x'; coef_y'];
end

%% Coefficients for pressures
yArray = -coef_y .* lengths;
d_yArray = [yArray(2:end); yArray(1)] - yArray;
xArray = coef_x .* lengths;
d_xArray = [xArray(2:end); xArray(1)] - xArray;

result2 = [d_yArray'; d_xArray'] .* 0.5;
end
