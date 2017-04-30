function [GCxy GCyx GCxy_bt GCyx_bt]=copu_gc_callfunc(data,mw,m_mir,h1,nlag_s,nlag_r,bt,nbt)

% Estimate the copula-based Granger causality between two time series
% 
% Input:
%       data    - ch x time
%       mw      - parameter of mirror processing to optimize the estimation 
%       m_mir   - grid size of copula
%       h1      - kernel bandwidth
%       nlag_s  - time lag of trigger signal
%       nlag_r  - time lag of receiver signal
%       bt      - flag for bootstrap
%       nbt     - number of bootstrap
% 
% Output:
%       GCxy    - Granger causality from 1st row to 2nd row in data
%       GCyx    - Granger causality from 2nd row to 1st row in data
%       GCxy_bt - boostrap Granger causality from 1st row to 2nd row in data
%       GCyx_bt - boostrap Granger causality from 2nd row to 1st row in data

% by Meng Hu @ Liang's lab at Drexel University

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

x=data(1,:);
y=data(2,:);
N=length(x);

%%%%% Source: x
%%%%% Receiver: y

for w=1:nlag_r

for l=1:nlag_s
    
    len_s=1:N-l;
    len_r=1:N;
    len_r_hist=1:N-w;
    len=min([length(len_s),length(len_r),length(len_r_hist)]);
    
    if length(len_s)>len
        len_s=len_s(1+length(len_s)-len:end);
    end
    
    if length(len_r)>len
        len_r=len_r(1+length(len_r)-len:end);
    end    
    
    if length(len_r_hist)>len
        len_r_hist=len_r_hist(1+length(len_r_hist)-len:end);
    end
        
    sample=[];
    sample(1,:)=x(len_s);
    sample(2,:)=y(len_r);
    sample(3,:)=y(len_r_hist); %% history     
    lsample=length(sample);    
    %%%%%%%%%%%%%%%%%%% conditional margin    

    Xquery4marg=1/lsample:1/lsample:1;
    Condquery=[1/m_mir : 2*1/m_mir : 1];
    [x1margcond x2margcond] = mar_cond_fast(sample',m_mir,h1,Xquery4marg,Condquery);
    %%%%%%%%%%%%%%%%%%% conditional margin 

    %%%%%%%%%%%%%%%%%%% conditional copula  

    Xquery4copu=[0 1/m_mir:1/m_mir:1];
    cdf_cond = ec_cdf_cond_fast(sample',m_mir,h1,Xquery4copu,Condquery,x1margcond,x2margcond);
    
    for q=1:size(cdf_cond,1)
        ec_cdf_tmp=squeeze(cdf_cond(q,:,:));
        ec_pdf=berncopupdf_mirror(ec_cdf_tmp,m_mir,mw); %% Bernstein copula pdf
        GCxy_tmp(q)=test4gc(ec_pdf,m_mir);
              
    end
    
    GCxy(w,l)=mean(GCxy_tmp); %% Receiver (self lag) X Source (Causal lag)

    %%%%%%%%%%%%%%%%%%% conditional copula   
    
    if bt==true    
    %%%%%%%%%%%%%%%%%%% Significance test
    %%%%%%%%%%%%%%%%%%%
    Condquery_bt=bootrsp(Condquery,nbt)';
    for kk=1:nbt
        tmp=Condquery_bt(kk,:);
        tmpp=repmat(tmp,[length(Condquery) 1]);
        tmppp=repmat(Condquery,[length(tmp) 1])';
        tmpindx=tmpp-tmppp;
        [tmpa,tmpb]=find(tmpindx==0);
        
        cdf_cond = ec_cdf_cond_boot_fast(sample',m_mir,h1,Xquery4copu,Condquery_bt(kk,:),x1margcond,x2margcond,tmpa);
        for q=1:size(cdf_cond,1)
            ec_cdf_tmp=squeeze(cdf_cond(q,:,:));
            ec_pdf=berncopupdf_mirror(ec_cdf_tmp,m_mir,mw); %% Bernstein copula pdf
            GCxy_tmp(q)=test4gc(ec_pdf,m_mir);
        end
        GCxy_bt(w,l,kk)=mean(GCxy_tmp);

    end    
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Significance test end    
    end

end

end
%%%%% x->y end





%%%%% Source: y
%%%%% Receiver: x

for w=1:nlag_r

for l=1:nlag_s
    
    len_s=1:N-l;
    len_r=1:N;
    len_r_hist=1:N-w;
    len=min([length(len_s),length(len_r),length(len_r_hist)]);
    
    if length(len_s)>len
        len_s=len_s(1+length(len_s)-len:end);
    end
    
    if length(len_r)>len
        len_r=len_r(1+length(len_r)-len:end);
    end    
    
    if length(len_r_hist)>len
        len_r_hist=len_r_hist(1+length(len_r_hist)-len:end);
    end
        
    sample=[];
    sample(1,:)=y(len_s);
    sample(2,:)=x(len_r);
    sample(3,:)=x(len_r_hist); %% history     
    lsample=length(sample);    
    %%%%%%%%%%%%%%%%%%% conditional margin    

    Xquery4marg=1/lsample:1/lsample:1;
    Condquery=[1/m_mir : 2*1/m_mir : 1];
    [x1margcond x2margcond] = mar_cond_fast(sample',m_mir,h1,Xquery4marg,Condquery);
    %%%%%%%%%%%%%%%%%%% conditional margin 

    %%%%%%%%%%%%%%%%%%% conditional copula  

    Xquery4copu=[0 1/m_mir:1/m_mir:1];
    cdf_cond = ec_cdf_cond_fast(sample',m_mir,h1,Xquery4copu,Condquery,x1margcond,x2margcond);
    
    for q=1:size(cdf_cond,1)
        ec_cdf_tmp=squeeze(cdf_cond(q,:,:));
        ec_pdf=berncopupdf_mirror(ec_cdf_tmp,m_mir,mw); %% Bernstein copula pdf
        GCyx_tmp(q)=test4gc(ec_pdf,m_mir);
        
        
    end
    
    GCyx(w,l)=mean(GCyx_tmp);
    %%%%%%%%%%%%%%%%%%% conditional copula          
    
    if bt==true    
    %%%%%%%%%%%%%%%%%%% Significance test
    %%%%%%%%%%%%%%%%%%%
    Condquery_bt=bootrsp(Condquery,nbt)';
    for kk=1:nbt
        tmp=Condquery_bt(kk,:);
        tmpp=repmat(tmp,[length(Condquery) 1]);
        tmppp=repmat(Condquery,[length(tmp) 1])';
        tmpindx=tmpp-tmppp;
        [tmpa,tmpb]=find(tmpindx==0);
        
        cdf_cond = ec_cdf_cond_boot_fast(sample',m_mir,h1,Xquery4copu,Condquery_bt(kk,:),x1margcond,x2margcond,tmpa);
        for q=1:size(cdf_cond,1)
            ec_cdf_tmp=squeeze(cdf_cond(q,:,:));
            ec_pdf=berncopupdf_mirror(ec_cdf_tmp,m_mir,mw); %% Bernstein copula pdf
            GCyx_tmp(q)=test4gc(ec_pdf,m_mir);
        end
        GCyx_bt(w,l,kk)=mean(GCyx_tmp);

    end    
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Significance test end    
    end
    
end

end
%%%%% x->y end




end


