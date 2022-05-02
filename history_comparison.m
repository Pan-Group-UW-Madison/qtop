clear all;
close all;
clc;

load('btop85_history');
load('gbd_history');

hold on;
plot(compliance, 'r.-');
plot(complianceHistory, 'b.-');
plot([0, length(complianceHistory)], [compliance(end), compliance(end)], 'k--');