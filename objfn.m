function J = objfn(P,N)

M = P(6*N+7:7*N+7);

tf = P(end);
J = tf;%
% J = -M(N+1);

end