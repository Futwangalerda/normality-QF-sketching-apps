%Table 4, 5
s1 = readmatrix("iidgauss1.csv");
t2 = readmatrix("srht.csv");
t3 = readmatrix("cs.csv");
s4 = readmatrix("sse1.csv");
t5 = readmatrix("haar.csv");
grid_m = 200:200:1600;
s1_mean=mean(s1,1);
s1_var=var(s1,0,1);
t2_mean=mean(t2,1);
t2_var=var(t2,0,1);
t3_mean=mean(t3,1);
t3_var=var(t3,0,1);
s4_mean=mean(s4,1);
s4_var=var(s4,0,1);
t5_mean=mean(t5,1);
t5_var=var(t5,0,1);

% sum=[grid_m; t1_mean; t2_mean; t3_mean; t4_mean; t5_mean; t1_var; t2_var; t3_var; t4_var; t5_var]';
% table_sum = array2table(sum);
% table_sum.Properties.VariableNames = {'sketch_size', 'iid', 'srht', 'countsketch', 'sse', 'haar', 'iidvar', 'srhtvar', 'countsketchvar', 'ssevar', 'haarvar'};
% writetable(table_sum,'sum.csv');