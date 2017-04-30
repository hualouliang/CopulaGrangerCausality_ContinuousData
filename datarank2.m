function rankx2=datarank2(x)

%%% Turn data into its rank
% Input:
%       x      - input data (time X channel)
% 
% Ouput:
%       rankx2 - rank of x

%%% Meng Hu @ Liang's lab at Drexel University
%%% 07/02/2012

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

[N d]=size(x);

%%%%%%%%%%%%%%%%%%%%%%%%% Transform real data into rank data
[y,idx]=sort(x,1);
rankx=[];
for n=1:d
    rankx(idx(:,n),n)=1:N; 
end
%%%%%%%%%%%%%%%%%%%%%%%%% Transform real data into rank data (end)

%%%%%%%%%%%%%%%%%%%%%%%%% Correct above rank results for the points with same values
rankx2=[];
for i=1:d
    tmp=rankx(:,i);
    k=1;
    for j=1:length(tmp)-1
        tmp1=find(tmp==k);
        tmp2=find(tmp==k+1);
        if abs(x(tmp1(1),i) - x(tmp2,i)) < 1e-10
           tmp(find(tmp>k))=tmp(find(tmp>k))-1;
        else
            k=k+1;
        end
    end
    
    tmpnew=tmp;
    for q=1:max(tmp)
        tmp1=find(tmp==q);
        tmpnew(find(tmp>q))=tmpnew(find(tmp>q))+length(tmp1)-1; 

    end
    
    rankx2(:,i)=tmpnew;
end
%%%%%%%%%%%%%%%%%%%%%%%%% Correct above rank results for the points with same values STOP


end