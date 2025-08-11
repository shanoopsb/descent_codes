function Pn = Legendre_poly(n)

if n == 0
   Pn = 0;
elseif n == 1
   Pn = ones(size(x));
else
   Pn_1 = ones(size(x));
   Pn   = x;

   for k = 1 : n-1
        PnP1 = 1/(k+1)*x*Pn - k*Pn_1;
        Pn_1 = Pn;
        Pn   = PnP1;
   end


end