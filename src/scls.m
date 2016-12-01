%% this is a linear least square solution realization
function a = linear_least_square (S, y)
	% S:		input for endmember hyperdata metrix
	% y:		input for one hyperdata consists of many band
	% a:		output for abundance solution for one hyperdata
    % the fomular for linear least square with constrains of sum of a is 1 is :
    %        _         		     _
	%		|	  (S'*S)^1_*1_'   |                 (S'*S)^*1_
	%	a=	| I - --------------  |*(S'*S)^*S'*y + --------------
	%		|_    1_'(S'*S)^*1_  _|                 1_'(S'*S)^1_
	%
	%

	[row_S, col_S] = size(S);

	STS = S' * S;
	STS_inv = inv (STS);
	one_vector = ones (col_S, 1);

	a1 = STS_inv * one_vector * one_vector' ;
	b1 = one_vector' * STS_inv * one_vector;

	b2 = a1 / b1;
	a2 = eye (col_S);

	a3 = a2 - b2;
	b3 = STS_inv * S' * y;

	a4 = STS_inv * one_vector;
	b4 = one_vector' * STS_inv * one_vector;

	a5 = a3 * b3;
	b5 = a4 / b4;

	a = a5 + b5; % the final solution
end
