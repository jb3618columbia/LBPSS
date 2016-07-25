function E = IsingEnergy(sigma)

% ISINGENERGY(SIGMA) calculates the energy (number of edges with unequal 
%   ends) of the Ising state SIGMA.
%   VERY SLOW!
[row col] = size(sigma);

E = 0; 
for i=1:row*col
    c(1) = mod(i-1,row)+1;
    c(2) = ceil(i/row);
    E = E + ( sigma(i)~=sigma(c(1), mod(c(2),col)+1) ) + ...
        ( sigma(i)~=sigma(mod(c(1),row)+1, c(2)) );
end
