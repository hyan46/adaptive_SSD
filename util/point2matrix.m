function [ X ] = point2matrix(W,y)

X = W*0;
X(W) = y;

end

