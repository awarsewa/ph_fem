function f = lagrange_poly_deriv(points, order)
    % Wikipedia has the 1st and 2nd derivative... 
    % https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivatives
    poly = cell(1, length(points));
    k = length(points);
    switch order 
        case 1
            for j=1:k
                poly{j} = @(x) 0;
                for i=1:k
                    if i ~= j
                        inner_poly = @(x) (1/(points(j)-points(i)));
                        for m = 1:k                  
                            if m ~= i && m ~= j
                                inner_poly = @(x) inner_poly(x) * ((x - points(m))/(points(j)-points(m)));
                            end
                        end
                        poly{j} = @(x) poly{j}(x) + inner_poly(x);
                    end
                end
            end
        case 2
            for j=1:k
                poly{j} = @(x) 0;
                for i=1:k
                    if i ~= j
                        inner_poly = @(x) 0;
                        for m=1:k
                            if m ~= j && m ~= i
                                inner_poly_2 = @(x) 1/(points(j)-points(m));
                                for p = 1:k                  
                                    if p ~= i && p ~= j && p ~= m
                                        inner_poly_2 = @(x) inner_poly_2(x) * (x - points(p))/(points(j)-points(p));
                                    end
                                end
                                inner_poly = @(x) inner_poly(x) + inner_poly_2(x);
                            end
                        end
                        poly{j} = @(x) poly{j}(x) + (1/(points(j)-points(i))) * inner_poly(x);
                    end
                end
            end
        otherwise
            % Sorry for being to lazy to implement higher orders
            error('order > 2 is not supported yet');
    end
    f = @(x) arrayfun(@(a,b) a{1, 1}(b),repmat(poly, length(x), 1), ones(length(x), length(poly)).*x);
end