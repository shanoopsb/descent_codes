function R = eci2tcr(phi,lam)

% phi latitude
% lam longitude


R = [-sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi); -sin(lam) cos(lam) 0 ; -cos(phi)*cos(lam) -cos(phi)*sin(lam) -sin(phi)];

end