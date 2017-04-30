function [x1mar_cond x2mar_cond]=mar_cond_fast(x,m,h,Xquery,Condquery)

% Calculate conditional mariginal CDF
% 
% Input:
%       x          - input data (time X channel)
%       m          - grid size of copula
%       h          - kernel bandwidth
%       Xquery     - Query samples for forming copula
%       Condquery  - Query samples for conditional variable
%       
% Output:
%       x1mar_cond - distribution of first column in data conditional on third column 
%       x2mar_cond - distribution of second column in data conditional on third column 

%%% Meng Hu @ Liang's lab at Drexel University
%%% 06/07/2012

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

[N d]=size(x);
Ncond=length(Condquery);

rankx2=datarank2(x);

Fx=rankx2/N;

h=h*1/m;

%%%%%%%%%%%%%%%%%% loop condition
% tic
x1mar_cond=[];
for p=1:Ncond   %% conditional query   
    condQueryRepeated = repmat(Condquery(p),N,1);
    Kcond = exp( -((condQueryRepeated - Fx(:,3)).^2) / h); %% Gaussian
    Kcondtmp=repmat(Kcond,1,N);
    lr=bsxfun(@le,Fx(:,1),Xquery);
    x1mar_cond(p,:)=sum(Kcondtmp.*lr,1)./sum(Kcondtmp,1); 
end

x2mar_cond=[];
for p=1:Ncond   %% conditional query   
    condQueryRepeated = repmat(Condquery(p),N,1);
    Kcond = exp( -((condQueryRepeated - Fx(:,3)).^2) / h); %% Gaussian
    Kcondtmp=repmat(Kcond,1,N);
    lr=bsxfun(@le,Fx(:,2),Xquery);
    x2mar_cond(p,:)=sum(Kcondtmp.*lr,1)./sum(Kcondtmp,1); 
end
% toc
%%%%%%%%%%%%%%%%%% loop condition


end
