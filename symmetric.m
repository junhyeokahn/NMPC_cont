function A = symmetric(X)

A = (X + X') - diag(diag(X));

end