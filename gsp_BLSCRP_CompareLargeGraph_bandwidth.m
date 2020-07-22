%    This code is for the  Fig.6(b) and Fig.7(b) in the paper 
%    Author:   Qian Zhang and Lihua Yang
%    Date:     1 July 2020
%    Version:  2020.7.01

%    varing the size of bandwidth, from n=1000 to n=7000
%    fixed N=20000, l=200

close all;
clear all;

G=gsp_sensor(20000);

nsize=[1000:1000:7000];

iNum=length(nsize); 

TIME=zeros(iNum,6);   %  6 lines

param.order = 100;
param.filter='lp-jch';
param.lk_est_method='fast';
param.gamma=0.1;


ell=200;


for i=1:iNum

n = nsize(i)

t1=clock;
[U1,D1]=eigs(G.L,n,'sm');
[D1,index] = sort(diag(D1),'ascend'); 
U1=U1(:,index);
t2=clock;

TIME(i,1)=etime(t2,t1);

t1=clock;
[Vn,D] = gsp_BLSCRP_Alg2(G,n,param);
t2=clock;
TIME(i,2)=etime(t2,t1);

loops=100;
Anal_time=0;
Alg1_time=0;

for j=1:loops
cutoffcoeff=rand(n,1);
cutoffcoeff=sort(cutoffcoeff,'descend');
f=U1*cutoffcoeff;             

p = randperm (G.N);          
x0 = p(1:ell);
y0 = f(x0);
U=U1(:,2:n);
mu=D1(2:n);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=clock;
y1=gsp_BLSCRP_Analytic(U1,D1,x0,y0,param);
t2=clock;
Anal_time=Anal_time+etime(t2,t1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1=clock;
y2=gsp_BLSCRP_Alg1(U,mu,x0,y0,param);
t2=clock;
Alg1_time=Alg1_time+etime(t2,t1);
end 
TIME(i,3)=Anal_time/loops;
TIME(i,4)=Alg1_time/loops;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TIME(i,5)=TIME(i,1)+TIME(i,3);
TIME(i,6)=TIME(i,2)+TIME(i,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

end

%%%%%%%%%%%%  Fig.1
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
plot(1:iNum,TIME(:,3),'b--o','LineWidth',2);
hold on;
plot(1:iNum,TIME(:,4),'r--o','LineWidth',2);
hold on;
 
set(gca,'xticklabel',{'1','2','3','4','5','6','7'},'linewidth',1.5);
xlabel('k','fontsize',18)
ylabel('TIMEs','fontsize',18)
set(gca,'FontSize',18);
legend( 'Anal','Alg1','Location','NorthWest'); 
%saveas(gcf, 'CompareLargeGraph_bandwidth1.jpg');
%%%%%%%%%%%%  Fig.2

figure('color',[1 1 1]);
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
 
plot(1:iNum,TIME(:,5),'b--o','LineWidth',2);
hold on;
plot(1:iNum,TIME(:,6),'r--o','LineWidth',2);
hold on;
 
set(gca,'xticklabel',{'1','2','3','4','5','6','7'},'linewidth',1.5);
xlabel('k','fontsize',18)
ylabel('TIMEs','fontsize',18)
set(gca,'FontSize',18);
legend( 'eigs+Alg1','Alg2+Alg1','Location','NorthWest'); 
%saveas(gcf, 'CompareLargeGraph_bandwidth2.jpg');

%save('data_CompareLargeGraph_bandwidth')
