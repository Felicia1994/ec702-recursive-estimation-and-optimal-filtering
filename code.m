% jingma
% 04/23/2018

clear;
[Y,T,w,b,profit] = initialize();
%plot_pro(w,b,profit);

A = 1;
%A = exp(w/252);
B = 0;
G = 1;
C = 1;

P_1 = 0; % P^(-1)
Z = 0;
X = zeros(T,1);
P = zeros(T,1);
temp = Y(1);
q = 0.0001;
r = 0.03;

for t = 1:T
    P_1 = P_1 + C * (r*temp)^(-2) * C;
    Z = Z + C * (r*temp)^(-2) * Y(t);
    M = A^(-1) * P_1 * A^(-1);
    N = M * G * ((q*temp)^(-2) + G * M * G)^(-1);
    P_1 = M - N * G * M;
    Z = (1 + M * G * (q*temp)^2 * G)^(-1) * A^(-1) * Z;
    temp = P_1^(-1) * Z;
    X(t) = temp;
    P(t) = P_1^(-1);
end

%avg_100 = cal_avg(Y,100);

figure;
plot(1:T,Y,'k','LineWidth',0.5);
xlabel('day','FontSize',14);
ylabel('daily stock price','FontSize',14);
hold on;
plot(1:T,X,'r','LineWidth',0.5);
legend({'real','filtered'},'Location','northwest');

figure;
fill([(1:T)';flipud((1:T)')],[X-sqrt(P);flipud(X+sqrt(P))],[.9 .9 .9],'linestyle','none');
xlabel('day','FontSize',14);
ylabel('filtered daily stock price','FontSize',14);
hold on;
plot(1:T,X,'r','LineWidth',0.5);

figure;
fill([(1:T)';flipud((1:T)')],[-sqrt(P)./X;flipud(sqrt(P)./X)],[.9 .9 .9],'linestyle','none');
xlabel('day','FontSize',14);
ylabel('normalized daily stock price residual','FontSize',14);
hold on;
plot(1:T,Y./X-1,'k','LineWidth',0.5);

%plot(1:T,X-Y,'r','LineWidth',1);
%plot(1:T,avg_100,'k','LineWidth',1);

%%%%%%%%%% functions %%%%%%%%%%
function [Y,T,w,b,profit] = initialize()

AMZN = readtable('AMZN.csv');
Y = flipud(AMZN(:,2));
Y = str2double(table2array(Y));
T = length(Y);
% revenue = [74452,88988,107006,135987,177866];
profit = [11686,15470,21945,30103,40683];

x = 1:5;
y = log(profit);
mu_x = mean(x);
mu_y = mean(y);
x = x - mu_x;
y = y - mu_y;
w = (x*y')/(x*x');
b = mu_y - w*mu_x;

end

function plot_pro(w,b,profit)

figure;
plot(2013:2017,profit,'-o','LineWidth',2);
xticks(2013:2017);
xlabel('year','FontSize',14);
ylabel('profit (million dollars)','FontSize',14);

figure;
plot(2013:2017,log(profit),'-o','LineWidth',2);
hold on;
xticks(2013:2017);
xlabel('year','FontSize',14);
ylabel('log of profit (million dollars)','FontSize',14);
plot(2013:2017,w*(1:5)+b,'-r','LineWidth',1);
legend({'real','fitting'},'Location','best');

end

function avg = cal_avg(x,k)

avg = zeros(length(x),1);

for i = 1:length(x)
    if i<k
        avg(i) = mean(x(1:i));
    else
        avg(i) = mean(x(i-k+1:i));
    end
end

end

function avg_exp = cal_avg_exp(x,rho)

avg_exp = zeros(length(x),1);
avg_exp(1) = x(1);

for i = 2:length(x)
    avg_exp(i) = (1-rho)*avg_exp(i-1) + rho*x(i);
end

end
