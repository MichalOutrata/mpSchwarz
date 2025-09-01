function val = b2( x,y, AdvectStrngth )

%val = AdvectStrngth.*(-1.*y.*(y-1).*(1-2.*x));
% val = -(x-0.5);
val = AdvectStrngth.*(1.*x.*(x-1).*(1-2.*y));

end

