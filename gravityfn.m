function [gx,gy,gz] = gravityfn(x,y,z,N,param)

rT = param.rT;
Mu = param.Mu;


gx = zeros(N+1,1);
gy = zeros(N+1,1);
gz = zeros(N+1,1);

for j = 1:N+1    
    r = [x(j);y(j);z(j)];
    rvec = r + rT;
    rnorm = norm(rvec);

    g     = -Mu*rvec/rnorm^3;

    gx(j,1) = g(1);
    gy(j,1) = g(2);
    gz(j,1) = g(3);

end



end