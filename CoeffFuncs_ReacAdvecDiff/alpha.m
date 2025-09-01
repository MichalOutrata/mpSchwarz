function val = alpha( x,y )

% if (y > 0.25) && (y < 0.50) && (x > 0.25) && (x < 0.50)
%   val = 1e+5;
% else
%   val = 1;
% end

if norm( [x-0.5,y-0.1] ) < 0.25
  val = 1e+6;
else
  val = 1;
end

% val = (x+y).^2.*exp(x-y);


end

