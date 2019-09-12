clc
clear all

A = csvread('results.csv');

plot(A(:,1),A(:,2))
xlabel('iterations')
ylabel('cost')