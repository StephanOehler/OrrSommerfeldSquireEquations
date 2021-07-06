function [ff, C] = bary_interp_new(x , xx , c , f)

N = length(x)-1;
m = length(xx);

% get the interpolating matrix C
C = zeros(m,N+1);
for j = 1:m % for each interpolant grid point (each row of C)
  % Barycentric formula as matrix
  diff = xx(j)-x;
  temp = c./diff;
  C(j,:) = temp/sum(temp);
  % deal with case where z(j)=(one of Chebyshev points)
  if min(abs(diff))==0
      C(j,:) = 0;      
      C(j,diff==0) = 1;
  end
end
ff = C * f;