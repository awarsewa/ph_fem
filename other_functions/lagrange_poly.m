function f = lagrange_poly(points)
    poly = cell(1, length(points));
    k = length(points);
    for j=1:k
        poly{j} = @(x) 1;
        for m=1:k
            if m ~= j
                poly{j} = @(x) poly{j}(x) * ((x - points(m))/(points(j)-points(m)));
            end
        end
    end
    f = @(x) arrayfun(@(a,b) a{1, 1}(b),repmat(poly, length(x), 1), ones(length(x), length(poly)).*x);
end