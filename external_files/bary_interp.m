function ff = bary_interp(x , xx , c , f)

N = length(x)-1;

numer = zeros(size(xx));
denom = zeros(size(xx));

exactx  = zeros(size(x ));
exactxx = zeros(size(xx));

for j = 1:N+1 % for each Chebyshev grid point
  xdiff = xx-x(j);
  temp = c(j)./xdiff;
  numer = numer + temp*f(j);
  denom = denom + temp;
  if min(abs(xdiff))==0, exactx(j) = 1; end
  exactxx(xdiff==0) = 1;
end
ff = numer./denom;

indx = find(exactx); indxx = find(exactxx);
ff(indxx) = f(indx);

%ff(find(exact)) = f(exact(find(exact)));
