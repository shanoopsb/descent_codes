function [x,D] = Dleg(N)


Np = N+1;
xc = cos(pi*(0:N)/N)';
xu = linspace(-1,1,Np)';

if N<3
    x = xc;
else
    x = sin(pi*xu)./(4*N) + xc;
end

xold = 2;

P = zeros(Np,Np);


while (max(abs(x-xold))> eps)
    
    xold = x;

    P(:,1) = 1;
    P(:,2) = x;

    for k = 2:N
        
        P(:,k+1) = 1/k*((2*k-1)*x.*P(:,k) - (k-1)*P(:,k-1));

    end

    x = xold - (x.*P(:,Np) - P(:,N))./(Np*P(:,Np));

end

X     = repmat(x,1,Np);
Xdiff = X - X' + eye(Np);

L                    = repmat(P(:,Np),1,Np);
L(1:Np+1:Np*Np)      = 1;

D               = L./(Xdiff.*L');
D(1:Np+1:Np*Np) = 0;
D(1,1)          =  N*Np/4;
D(Np,Np)        = -N*Np/4;

end