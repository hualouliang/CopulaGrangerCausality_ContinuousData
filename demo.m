% An example to demonstrate how to use our code to perform the Granger causality

% Please cite the following paper if you use this software:
% "Hu & Liang, A copula approach to assessing Granger causality, NeuroImage, 2014."

% This example can be deemed as a special case (q = 0.5) in simulation Model V in above reference.

% Meng Hu @ Liang's lab at Drexel University

clear

dl=1000; %% data length
c=0.5; %% effect coef of causality

%% data generation (causal relationship of nonlinear plus variance)

% X->Y
x=[];
y=[];
x(1)=randn(1);
y(1)=randn(1);
for n=1:dl-1
    x(n+1)=0.2*x(n)+randn(1);
    y(n+1)=0.1*y(n)+c*cos(x(n))*exp((-(x(n)).^2-(y(n)).^2)/8)+sqrt(0.2*y(n).^2+(1-c)*x(n).^2)*randn(1);
end      

%% Copula GC

data=[];
data(1,:)=x;
data(2,:)=y;

mw=2;
m_mir=20;
h1=3;
nlag_s=1;
nlag_r=1;
bt=true;
nbt=50;

[GCxy GCyx GCxy_bt GCyx_bt]=copu_gc_callfunc(data,mw,m_mir,h1,nlag_s,nlag_r,bt,nbt);

alpha=0.05;
GCxy_bt=sort(squeeze(GCxy_bt));
GCyx_bt=sort(squeeze(GCyx_bt));

if GCxy > GCxy_bt(fix(nbt*(1-alpha)))
    fprintf('True positive \n'); % correctly identify GC
end
if GCyx > GCyx_bt(fix(nbt*(1-alpha)))
    fprintf('False positive \n'); % falsely identify GC
end
