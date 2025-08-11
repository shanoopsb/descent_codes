function xout = Linterp(x_values,y_values,x)

xout = 0;
n = length(x_values);

for i = 1:n
    l = 1;
    for j = 1:n

        if i ~= j
           l = l.*((x - x_values(j))./(x_values(i)-x_values(j)));
        end
    end
xout = xout + y_values(i).*l;
end

end