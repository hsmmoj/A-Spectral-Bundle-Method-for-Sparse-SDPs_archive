function [vecA, vecB, vecC] = Vectorization(A,B,C, )

vecA = reshape(A, [size(A,1)*size(A,2), size(A,3)])';
vecB = B;
vecC = reahpe(C, [size(C,1)*size(C,2), size(C,3)])';


end