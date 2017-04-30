function [cdf_cond m]=ec_cdf_cond_boot_fast(x,m,h,Xquery,Condquery,x1margcond,x2margcond,tmpa)

% Estimate the bootstrap of cdf of conditional copula
% Please see ec_cdf_cond_fast.m

%%% Meng Hu @ Liang's lab at Drexel University

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

[N d]=size(x);
Nxq=length(Xquery);

rankx2=datarank2(x);

Fx=rankx2/N;

h=h*1/m;

cdf_cond=[];
for p=1:length(Condquery)   %% conditional query
    
    x1marg=x1margcond(tmpa(p),:)';
    x2marg=x2margcond(tmpa(p),:)';
    
    Fcondx1=x1marg(randperm(N)); %% randomly sample
    Fcondx2=x2marg(randperm(N));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% no loop for cdf
    condQueryRepeated = repmat(Condquery(p),N,1);
    Kcond = exp( -((condQueryRepeated - Fx(:,3)).^2) / h); %% Gaussian
    Kcondtmp=repmat(Kcond,1,Nxq.^2);
    
    xqtmp=meshgrid(Xquery);
    xqtmp1=reshape(xqtmp,1,Nxq.^2);
    xqtmp2=repmat(Xquery,1,Nxq);
    
    lr1=bsxfun(@le,Fcondx1,xqtmp1);
    lr2=bsxfun(@le,Fcondx2,xqtmp2);
    
    cdftmp=sum(Kcondtmp.*lr1.*lr2,1)./sum(Kcondtmp,1);
    cdf_cond(p,:,:)=reshape(cdftmp,Nxq,Nxq)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% no loop for cdf
    
end
    
end
