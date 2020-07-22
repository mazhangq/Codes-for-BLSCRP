% experiment on RS-movielens100k

% the datasets "ml-100k" can be found in "https://grouplens.org/"

clc;
clear all;
close all;

load('D:\ml-100k\mat\u1base.mat');
u{1}.base=u1base;
clear u1base;
load('D:\ml-100k\mat\u1test.mat');
u{1}.test=u1test;
clear u1test;
load('D:\ml-100k\mat\u2base.mat');
u{2}.base=u2base;
clear u2base;
load('D:\ml-100k\mat\u2test.mat');
u{2}.test=u2test;
clear u2test;
load('D:\ml-100k\mat\u3base.mat');
u{3}.base=u3base;
clear u3base;
load('D:\ml-100k\mat\u3test.mat');
u{3}.test=u3test;
clear u3test;
load('D:\ml-100k\mat\u4base.mat');
u{4}.base=u4base;
clear u4base;
load('D:\ml-100k\mat\u4test.mat');
u{4}.test=u4test;
clear u4test;
load('D:\ml-100k\mat\u5base.mat');
u{5}.base=u5base;
clear u5base;
load('D:\ml-100k\mat\u5test.mat');
u{5}.test=u5test;
clear u5test;
load('D:\ml-100k\mat\uabase.mat');
u{6}.base=uabase;
clear uabase;
load('D:\ml-100k\mat\uatest.mat');
u{6}.test=uatest;
clear uatest;
load('D:\ml-100k\mat\ubbase.mat');
u{7}.base=ubbase;
clear ubbase;
load('D:\ml-100k\mat\ubtest.mat');
u{7}.test=ubtest;
clear ubtest;

U=943;
I=1682;
 

for r=1:5
sim=1;
ubase=u{r}.base;
utest=u{r}.test;
M1=zeros(U,I);
M1((ubase(:,2)-1)*U+ubase(:,1))=ubase(:,3);    
testM1=zeros(U,I);
testM1((utest(:,2)-1)*U+utest(:,1))=utest(:,3);    

M=M1;
testM=testM1;
mask=(M>0);

 

labelM=[];
labeltestM=[];
for i=1:I
    labelM=[labelM,length(find(mask(:,i)))];
    %labeltestM=[labeltestM,length(find(testmask(:,i)))];
end 
II=find(labelM>30);
IIlength=length(II);

M=[];
newM=[];
testM=[];
mask=[];
M=M1(:,II);
testM=testM1(:,II);
mask=(M>0);
testmask=(testM>0);
numtest=sum(sum(testM~=0));   

%[aM,mU,mUI]=AdjustUI(M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errs=[];
if sim==3
    Wi=Simxy(aM,mask,1);
else
    Wi=Simxy(M,mask,sim);
end
D=diag(sum(Wi,2));
L=D-Wi;
G=gsp_graph(Wi);

N = G.N; 
W=  G.W;
beta = 1;  
param.gamma = 0.1;  
n=round(N/3);  


t1=clock;
 
[eigU,eigD]=eig(L);
eigD=diag(eigD);
[eigD,index] = sort(eigD,'ascend'); 
eigU=eigU(:,index);

mu=eigD(2:n);
Un=eigU(:,2:n);
 
 

t2=clock;
eigtime=etime(t2,t1)

D1=eigD(1:n);
U1=eigU(:,1:n);

t1=clock;

for j=1:IIlength
    
S=find(mask(:,j));

f=M(:,j);

y=f(S);
 
x0=S;

y0=y;

ell=numel(x0);
 
B=Un(x0,:);
 
temp=repmat(1./mu,1,size(B,1));

temp=temp.*B';

Tl=B*temp;
 
temp2=Tl*ones(ell,1);

GG=param.gamma*Tl+Tl*Tl-1/ell*(temp2*temp2');

d=Tl*y0-1/ell*temp2*(ones(1,ell)*y0);

xi=linsolve(GG,d);

g=Un*(temp*xi);
 
temp3=- sum(g(x0)-y0)/ell;

ff=temp3+g;
 
newM(:,j)=ff;

end

t2=clock;
fillingtime=etime(t2,t1)
totaltime=fillingtime+eigtime

M2=newM.*testmask;
A=abs(M2-testM);
error=0;

for j=1:IIlength
    error=error+norm(A(:,j),1);
end 
error=error/numtest
end

 