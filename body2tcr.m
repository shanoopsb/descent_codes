function R = body2tcr(phi,theta,psi)


R(1,1) =  cos(theta)*cos(psi);
R(2,1) =  cos(theta)*sin(psi);
R(3,1) = -sin(theta);

R(1,2) =  sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
R(2,2) =  sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi);
R(3,2) =  sin(phi)*cos(theta);

R(1,3) =  cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
R(2,3) =  cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
R(3,3) =  cos(phi)*cos(theta);

end