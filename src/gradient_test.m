Y = [2; 5];
C = [1,1; 2, 3];
X_INIT = [0.5; 0.5];
RATE = 0.01;
TOLERANT = 0.0001;

[X, E] = linear_least_square_gradient (Y', C', X_INIT', RATE, TOLERANT)