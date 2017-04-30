function gc=test4gc(c,nbin)

% Summation of copula
% Input:
%      c    - estimated copula
%      nbin - grid size of copula
%      
% Output:
%      gc   - summation of copula

%%%% Meng Hu @ Liang's lab at Drexel University

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

x=c(:);
xind=find(x>0);
gc=sum(log(1./(x(xind))))/nbin.^2;

end