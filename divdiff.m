function TDD = divdiff(X, Y)
%   Newton's Method for Divided Differences.
%   The following formula is solved:
%       Pn(x) = f(x0) + f[x0,x1](x-x0) + f[x0,x1,x2](x-x0)(x-x1) + ...
%           + f[x0,x1,..,xn](x-x0)(x-x1)..(x-x[n-1])
%       where f[x0,x1] = (f(x1-f(x0))/(x1-x0)
%             f[x0,x1,..,xn] = (f[x1,..,xn]-f[x0,..,x_[n-1]])/(xn-x0)
%	X = [ x0 x1 .. xn ]
%	Y = [ y0 y1 .. yn ] 
    if nargin ~= 2 %check the validity of the function parameters
        error('divdiff: invalid input parameters'); 
    end
    [ p, m ] = size(X); % m points, polynomial order <= m-1
    if p ~= 1 || p ~=size(Y, 1) || m ~= size(Y, 2) %check that both vectors have equal dimensions
        error('divdiff: input vectors must have the same dimension'); 
    end
    TDD = zeros(m, m); %create a zero matrix where the difference values will be added to it 
    TDD(:, 1) = Y'; %enter values of y as the first column
    for j = 2 : m %loop over the size of the table to get the y_coordinates
        for i = 1 : (m - j + 1)  %loop for the x coorddinates
            TDD(i,j) = (TDD(i + 1, j - 1) - TDD(i, j - 1)) / (X(i + j - 1) - X(i)); %recursively finding the divided table
        end
    end
end