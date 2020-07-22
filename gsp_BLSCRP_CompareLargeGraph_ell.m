%    This code is for the  Fig.6(c)  in the paper 
%    Author:   Qian Zhang and Lihua Yang
%    Date:      1 July 2020
%    Version:  2020.7.01 

%   varing the size of known labeled points, from l=300k-100, k=1,...,7.
%   fixed N=20000, n=7000

close all;
clear all;

load('data_CompareLargeGraph_bandwidth')
% 
ellsize=[200:300:2000];

iNum=length(ellsize); 

TIME=zeros(iNum,6);   %  6 lines

param.order = 100;
param.filter='lp-jch';
param.lk_est_method='fast';
param.gamma=0.1;
p = randperm (G.N);          

for i=1:iNum

ell = ellsize(i)

x0 = p(1:ell);
y0 = f(x0);
U=U1(:,2:n);
mu=D1(2:n);

% t1=clock;
% [U1,D1]=eigs(G.L,n,'sm');
% [D1,index] = sort(diag(D1),'ascend'); 
% U1=U1(:,index);
% t2=clock;
% 
% TIME(i,1)=etime(t2,t1);
% 
% 
%  
% t1=clock;
% [Vn,D] = gsp_BLSCRP_Alg2(G,n,param);
% t2=clock;
% TIME(i,2)=etime(t2,t1);

loops=100;
Anal_time=0;
Alg1_time=0;

for j=1:loops
cutoffcoeff=rand(n,1);
% cutoffcoeff=sort(cutoffcoeff,'descend');
f=U1*cutoffcoeff;            

 
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
% TIME(i,5)=TIME(i,1)+TIME(i,3);
% TIME(i,6)=TIME(i,2)+TIME(i,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

end


load('data_CompareLargeGraph_ell')
figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
plot(1:iNum,TIME(:,3),'b--o','LineWidth',2);
hold on;
plot(1:iNum,TIME(:,4),'r--o','LineWidth',2);
hold on;
axis([-inf,inf,0,16])
set(gca,'xticklabel',{'1','2','3','4','5','6','7'},'linewidth',1.5);
xlabel('k','fontsize',18)
ylabel('TIMEs','fontsize',18)
set(gca,'FontSize',18);
legend( 'Anal','Alg1','Location','NorthWest'); 
saveas(gcf, 'CompareLargeGraph_ell.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TIME(:,3)=fliplr(TIME(:,3));
% TIME(:,4)=fliplr(TIME(:,4));

t1=sort(TIME(:,3),'descend');
t2=sort(TIME(:,4),'descend');

figure
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
plot(1:iNum,t1/t1(1),'b--o','LineWidth',2);
hold on;
plot(1:iNum,t2/t2(1),'r--o','LineWidth',2);
hold on;
 
set(gca,'xticklabel',{'1','2','3','4','5','6','7'},'linewidth',1.5);
xlabel('k','fontsize',18)
ylabel('TIMEs','fontsize',18)
set(gca,'FontSize',18);
legend( 'Anal','Alg1','Location','NorthWest'); 
saveas(gcf, 'CompareLargeGraph_ellratio.jpg');




%save('data_CompareLargeGraph_ell')