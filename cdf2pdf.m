function pdf=cdf2pdf(ec_cdf)

% Convert cdf to pdf
% Input:
%       ec_cdf    - input CDF (square matrix)
%       
% Output:
%       pdf       - estimated PDF

%%% Meng Hu @ Liang's lab at Drexel University

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."


[n1,n2]=size(ec_cdf);
        
ec_cdf00=ec_cdf(1:n1-1,1:n1-1);
            
ec_cdf10=ec_cdf(1:n1-1,2:n1);
               
ec_cdf01=ec_cdf(2:n1,1:n1-1);
           
ec_cdf11=ec_cdf(2:n1,2:n1);
          
pdf = ec_cdf00 + ec_cdf11 - ec_cdf10 - ec_cdf01;

end