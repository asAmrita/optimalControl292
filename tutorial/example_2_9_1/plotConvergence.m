clc
clear all

%%
A = csvread('error.csv');
plot(A([1:end],1),A([1:end],2))
xlabel('iterations')
ylabel('||u-u_d||_2')


%%
A = csvread('cost.csv');
plot(A([2:end],1),A([2:end],2))
xlabel('iterations')
ylabel('cost')