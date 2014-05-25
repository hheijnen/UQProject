%Write here your likelihood!
function [fval,output]=uq14_model(x,data)


%Return SSE (not important)
output=-0.5*(sum(x(1:(end)).^2)); 



%mu=[0,0];
%si=[1 0;0 1];
%mu2=[-3,-3];
%si2=[0.1 0;0 1];


%Return the value of the LOG-LIKELIHOOD- VERY IMPORTANT!!!
%fval = log(mvnpdf(x,mu,si)+ mvnpdf(x,mu2,si2));
fval = SIRss(x,data.expvar);















%%%%%%%%%%Function to get the Log-Determinant of a Matrix%%%%%%%%%%%%%%%%%%
function v = logdet(A, op)

assert(isfloat(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
    'logdet:invalidarg', ...
    'A should be a square matrix of double or single class.');

if nargin < 2
    use_chol = 0;
else
    assert(strcmpi(op, 'chol'), ...
        'logdet:invalidarg', ...
        'The second argument can only be a string ''chol'' if it is specified.');
    use_chol = 1;
end

if use_chol
    v = 2 * sum(log(diag(chol(A))));
else
    [L, U, P] = lu(A);
    du = diag(U);
    c = det(P) * prod(sign(du));
    v = log(c) + sum(log(abs(du)));
end

