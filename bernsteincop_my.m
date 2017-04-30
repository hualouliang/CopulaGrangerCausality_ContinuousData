function b = bernsteincop_my(c,g)

% Bernstein approximation of copula
% 
% Input:
%      c    - copula
%      g    - grid size of copula
% 
% Output:
%       b    - Bernstein approximation

% Modified by Meng Hu @ Liang's lab at Drexel University, based on the code
% written by Robert Kopocinski, Wroclaw University of Technology,
%   for Master Thesis: "Simulating dependent random variables using copulas.
%   Applications to Finance and Insurance".

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

[m n] = size(c);

b=zeros(g);

for i=1:g
    for j=1:g
        b(i,j)=sum(sum( bernstein(i/g,1:m,m)'*(bernstein(j/g,1:m,m)).*c ));
    end
%     i
end

%------------------------------------------------------------
function y = bernstein(u,i,n)
% BERNSTEIN calculates values of Bernstein polynomials.
y=exp(gammaln(n+1)-gammaln(i+1)-gammaln(n-i+1)).*u.^i.*(1-u).^(n-i);

