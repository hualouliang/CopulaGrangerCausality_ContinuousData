function ec_pdf_smooth=berncopupdf_mirror(ec_cdf,m,mw)

% Estimate the copula pdf based on Bernstein
% 
% Input:
%      ec_cdf   - copula cdf
%      m        - grid size of copula
%      mw       - parameter of mirror processing to optimize the estimation
% Output:
%      copula pdf

%%% Meng Hu @ Liang's lab at Drexel University

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

delta_cdf = cdf2pdf(ec_cdf);
        
flip1=fliplr(delta_cdf);
flip2=flipud(delta_cdf);
flip3=flipud(flip1);
flip4=fliplr(flip2);

[u u]=size(delta_cdf);

bigpdf=[];

bigpdf(1:u,u+1:2*u)=flip2;
bigpdf(1+u:2*u,1:u)=flip1;
bigpdf(1+u:2*u,1+u:2*u)=delta_cdf;
bigpdf(1+u:2*u,1+2*u:3*u)=flip1;
bigpdf(1+2*u:3*u,u+1:2*u)=flip2;
bigpdf(1:u,1:u)=flip3;
bigpdf(1+2*u:3*u,1:u)=flip3;
bigpdf(1:u,1+2*u:3*u)=flip4;
bigpdf(1+2*u:3*u,1+2*u:3*u)=flip4;

pdf_mir=bigpdf(u-mw+1:2*u+mw,u-mw+1:2*u+mw);

pdf_mir_smooth = bernsteincop_my(pdf_mir,u+2*mw); %% Bernstein approximation

ec_pdf_smooth = m.^2*pdf_mir_smooth(mw+1:mw+u,mw+1:mw+u);


end