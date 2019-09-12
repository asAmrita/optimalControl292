clc
clear all

%%
A = csvread('convergence_example_harbir.csv');
loglog(A([1:end],1),A([1:end],2))

xlabel('num of cells')


%%
B = csvread('convergence_example_book_ex291.csv');
loglog(B([1:end],1),B([1:end],2))


xlabel('num of cells')