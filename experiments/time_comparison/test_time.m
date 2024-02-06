%Compare Computational Cost in Different Sketching

n = 2^11;
p = 15;
alpha = 0.1;
k = 1;
sim = 500;
zeta = 8;
grid_m = 200:200:1600;
c = [1, zeros(1, p-1)];
rng(15414514);
%Case 1
d = 1./(1:p);
D = diag(d);
O1 = randn(n, p);
[W, ~, ~] = svd(O1, 'econ');
O2 = randn(p, p);
[~, ~, V] = svd(randn(p, p), 'econ');
X = W * D * V';
%{
%Case 3
X1 = randn(n/2, p);
X2 = normrnd(5, 1, n/2, p);
X = [X1; X2];
%}

t1 = zeros(sim, length(grid_m));
r_gauss = zeros(sim, length(grid_m));
projvec_gauss = zeros(sim, length(grid_m));

for j = 1:length(grid_m)
    for i = 1:sim
        mt = floor(grid_m(j)^(1));
        S = randn(mt, n) / sqrt(mt); % S: m * n
        start_time_1 = datetime('now');
        [r_gauss(i, j), projvec_gauss(i, j)] = sketchingMethods.gauss(k, X, c, S);
        end_time_1 = datetime('now');
        t1(i, j) = vpa(seconds(end_time_1 - start_time_1), 4);
    end
end
writematrix(t1, 'iidgauss1.csv')

t2 = zeros(sim, length(grid_m));
r_srht = zeros(sim, length(grid_m));
projvec_srht = zeros(sim, length(grid_m));

for j = 1:length(grid_m)
    for i = 1:sim
        start_time_2 = datetime('now');
        [r_srht(i, j), projvec_srht(i, j)] = sketchingMethods.srht(k, grid_m(j), X, c);
        end_time_2 = datetime('now');
        t2(i, j) = vpa(seconds(end_time_2 - start_time_2), 4);
    end
end
writematrix(t2, 'srht.csv')

t3 = zeros(sim, length(grid_m));
r_cs = zeros(sim, length(grid_m));
projvec_cs = zeros(sim, length(grid_m));

for j = 1:length(grid_m)
    for i = 1:sim
        start_time_3 = datetime('now');
        [r_cs(i, j), projvec_cs(i, j)] = sketchingMethods.sse(k, grid_m(j), X, 1, c);
        end_time_3 = datetime('now');
        t3(i, j) = vpa(seconds(end_time_3 - start_time_3), 4);
    end
end
writematrix(t3, 'cs.csv')

t4 = zeros(sim, length(grid_m));
r_sse = zeros(sim, length(grid_m));
projvec_sse = zeros(sim, length(grid_m));

for j = 1:length(grid_m)
    for i = 1:sim
        start_time_4 = datetime('now');
        [r_sse(i, j), projvec_sse(i, j)] = sketchingMethods.sse(k, grid_m(j), X, zeta, c);
        end_time_4 = datetime('now');
        t4(i, j) = vpa(seconds(end_time_4 - start_time_4), 4);
    end
end
writematrix(t4, 'sse1.csv')

t5 = zeros(sim, length(grid_m));
r_haar = zeros(sim, length(grid_m));
projvec_haar = zeros(sim, length(grid_m));

for j = 1:length(grid_m)
    for i = 1:sim
        start_time_5 = datetime('now');
        [r_haar(i, j), projvec_haar(i, j)] = sketchingMethods.orth(k, grid_m(j), X, c);
        end_time_5 = datetime('now');
        t5(i, j) = vpa(seconds(end_time_5 - start_time_5), 4);
    end
end
writematrix(t5, 'haar.csv')
