

function result = projectionRatioAround(center,vertices)

projected = edgeLengthAround(center,vertices);

actual = sqrt(sum((vertices - center) .^ 2,2));

result = projected ./ actual;
end
