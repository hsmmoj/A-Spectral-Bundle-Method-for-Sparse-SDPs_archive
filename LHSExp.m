function bigMatrix = LHSExp(matrix, clique, dim)
m = size(matrix,2);
bigMatrix = zeros(dim,m);
bigMatrix(clique, :) = matrix;
end