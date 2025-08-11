% Example: LGR nodes (manually defined for illustration, replace with actual nodes)

N = 50;
x = Radau_roots(N);

D = compute_Lagrange_D(x);
disp(D);




function D = compute_Lagrange_D(x)
    % Computes the differentiation matrix D_{k,i} = L_i'(x_k)
    % where L_i(t) is the i-th Lagrange basis polynomial

    N = length(x) - 1;
    D = zeros(N+1, N+1);

    for k = 1:N+1
        for i = 1:N+1
            if k ~= i
                sum_val = 0;
                for l = 1:N+1
                    if l ~= i
                        prod_val = 1;
                        for j = 1:N+1
                            if j ~= i && j ~= l
                                prod_val = prod_val * (x(k) - x(j)) / (x(i) - x(j));
                            end
                        end
                        sum_val = sum_val + (1 / (x(i) - x(l))) * prod_val;
                    end
                end
                D(k,i) = sum_val;
            else
                % Diagonal term: L_i'(x_i)
                sum_val = 0;
                for j = 1:N+1
                    if j ~= i
                        sum_val = sum_val + 1 / (x(i) - x(j));
                    end
                end
                D(k,i) = sum_val;
            end
        end
    end
end
