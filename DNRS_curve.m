


%%
clear;
clc;
close all;
[hESCs_to_neurons_DNRS_matrix,genes] = xlsread('hESCs_to_neurons_DNRS_matrix.xlsx');
hESCs_to_neurons_DNRS_matrix_size=size(hESCs_to_neurons_DNRS_matrix);
for t=1:5
    [sort_result_curve,idx]=sort(hESCs_to_neurons_DNRS_matrix(:,t),'descend');
    tmp=sort_result_curve(1:196);
    case_result(t)=mean(tmp); 
end
plot(1:5,case_result,'r-*','LineWidth',3);
ylim([9 10.5])








