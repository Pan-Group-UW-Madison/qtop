close all;
clear all;
clc;

T1 = [0.40 0.86 2.30 6.96 8.83 6.51];
iter1 = [37 35 48 86 87 52];
t1 = T1 ./ iter1;

% T2 = [0.81 3.35 14.57 21.91 48.92];
% iter2 = [16 15 15 14 17];
% t2 = T2 ./ iter2;

N = [120*40 240*80 360*120 480*160 540*180 600*200];
hold on;
plot(log(N)/log(10), log(t1)/log(10), 'o-');
% plot(sqrt(N), log(t2)/log(10), 's-');

xlabel('$\log_{10}n_{\rho}$', 'Interpreter', 'latex');
ylabel('$\log_{10}t_{\mathrm{single}}$', 'Interpreter', 'latex');

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20)