% Analysis of the simulation results for the explosive rigidity percolation
% in mxn Miura-ori structures (the rectangular case)
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

% Pattern size (m*n)
mn_all = [25, 36, 49, 64, 81, 100, 121, 144, 169, 196, ...
    225, 256, 289, 324, 361, 400]; 

% number of choices in each step
k_all = [1, 2, 4, 8, 16, 32]; 

% selection rule (1: most efficient, 2: least efficient)
% rule_all = [1,2]; 
rule_all = [1]; 

% number of simulations for each setup
n_sim = 500;

%% Extract all results

% Initialize data storage
result_all = [];

for rule = 1 
    for k = k_all
        id_k = find(k_all == k);

        log_m_div_n_all = [];
        rho_w_all = [];
        
        for mn = mn_all
        
            num_all = 2:mn-1;
            % get all combinations of m*n with m,n>=2
            L_all = num_all(rem(mn,num_all)==0);
            
            m_all = zeros(length(L_all), 1);
            n_all = zeros(length(L_all), 1);

            for L = L_all
                
                id = find(L_all == L);

                load(['simulation_results/rec_', num2str(mn), ...
                    '/sim_rect_', num2str(mn), '_L_', num2str(L), '_k_', ...
                    num2str(k), '_rule_', num2str(rule), '.mat']);

                m = L;
                n = mn / L;
                M = m + 1;
                N = n + 1;
                [X, Y, Z] = generate_miura_ori(M, N);
                A_initial = RigidityMatrix(X, Y, Z, m, n);
                DoF_initial = 3 * M * N - calc_rank(A_initial) - 6;

                % percentage of planarity constraints added
                r = linspace(0, 1, m * n + 1)';

                % probability of getting a 1-DoF structure
                P = sum([DoF_initial * ones(n_sim, 1), dof_all] == 1, 1) / n_sim;

                % critical rho^*
                idx = find(P >= 0.5, 1); 
                r_star = r(idx);

                % rho^0 and rho_1
                idx_r0 = find(P == 0, 1, 'last'); 
                idx_r1 = find(P == 1, 1); 
                r0 = r(idx_r0);
                r1 = r(idx_r1);

                result_all = [result_all; m, n, k, r0, r1, r_star];
            end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyzing r_diff
m = result_all(:,1);
n = result_all(:,2);
k = result_all(:,3);
r_star = result_all(:,6);
r_min = (2*m+2*n-4)./(m.*n);
r_diff = r_star-r_min; % rho_diff

% plot the simulated rho_diff vs k for each L
f = @(k) (1+7.8*k./(sqrt(m.*n))).^(-1.4).*(1-(log(m./n)./log(m.*n/4)).^2);
fit_rdiff = f(k);

figure;
plot(fit_rdiff, r_diff,'o', 'Color', [201 0 22]/255, 'MarkerFaceColor', [201 0 22]/255,'LineWidth', 2);
xlabel('Fitted \rho_{diff}');
ylabel('Simulated \rho_{diff}');
hold on;
axis equal
plot([0 0.75], [0, 0.75], '-', 'Color', [201 0 22]/255,'LineWidth', 1);
axis([0 0.75 0 0.75])
set(gca, 'FontSize', 24);
set(gca, 'LineWidth', 2);
box on;

%% analyzing the results for a specific mn
m = result_all(:,1);
n = result_all(:,2);
k = result_all(:,3);
r_star = result_all(:,6);
r_min = (2*m+2*n-4)./(m.*n);
r_diff = r_star-r_min; % rho_diff

% plot the simulated rho_diff for a specific mn
cmap = color_scheme(length(k_all));
mn = 100;
figure;
hold on;
for iii = 1:length(k_all)
    id = find(m.*n == mn & k == k_all(iii));
    plot(log(m(id)./n(id))', r_diff(id)','-o','Color',cmap(iii,:),'MarkerFaceColor',cmap(iii,:),...
        'LineWidth', 2, 'DisplayName', ['k = ', num2str(k_all(iii))]);
end
legend show
xlabel('log(m/n)');
ylabel('Simulated \rho_{diff}')
title(['mn = ', num2str(mn)])
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([-5 5 0 1])
box on;

% plot the fitted rho_diff
k_range = 1:0.1:32;
mn = 100;
m = 2:0.1:mn/2;
n = mn./m;
figure;
hold on;
for iii = 1:length(k_all)
    kk = k_all(iii);
    f = (1+7.8*kk./(sqrt(m.*n))).^(-1.4).*(1-(log(m./n)./log(m.*n/4)).^2);
    plot(log(m./n)', f','-','Color',cmap(iii,:),'LineWidth', 2, 'DisplayName', ['k = ', num2str(kk)]);
end
legend show
xlabel('log(m/n)');
ylabel('Fitted \rho_{diff}')
title(['mn = ', num2str(mn)])
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([-5 5 0 1])
box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyzing rho_w
m = result_all(:,1);
n = result_all(:,2);
k = result_all(:,3);
r0 = result_all(:,4);
r1 = result_all(:,5);
r_min = (2*m+2*n-4)./(m.*n);
r_w = r1-r0; % rho_width

% fitting model for rho_w
f0 = @(k) (r_min-1./(m.*n)) + (1-(r_min-1./(m.*n))).*exp(-6.3*sqrt(min(k./sqrt(m.*n),sqrt(m.*n))));
f1 = @(k) (r_min)  + (1-r_min).*exp(-0.15*min(k,m.*n).^(1.2)./sqrt(sqrt(m.*n)));
f = @(k) (f1(k)-f0(k));

% Alternative formula proposed in SI
% L = sqrt(m.*n);
% rmin = (4*L-4)./L.^2;
% fit_r0 = (rmin-1./L.^2) + (1-(rmin-1./L.^2)).*exp(-6.3*sqrt(min(k./L,L)));
% fit_r1 = (rmin) + (1-rmin).*exp(-0.15*min(k,L.^2).^(1.2)./sqrt(L));
% f = @(k) (f1(k)-f0(k) - 1./(m.*n)).*(1-((log(m./n))./log(m.*n/4)).^2) + 1./(m.*n);

fit = f(k);

figure;
% plot fitted rho_w and simulated rho_w
plot(fit, r_w,'o', 'Color', [201 0 22]/255, 'MarkerFaceColor', [201 0 22]/255,'LineWidth', 2);
xlabel('Fitted \rho_w');
ylabel('Simulated \rho_w');
hold on;
% plot y = x
plot([0 0.75], [0, 0.75], '-', 'Color', [201 0 22]/255,'LineWidth', 1);
axis equal
axis([0 0.75 0 0.75])
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

%% analyzing the results for a specific mn
m = result_all(:,1);
n = result_all(:,2);
k = result_all(:,3);
r0 = result_all(:,4);
r1 = result_all(:,5);
r_w = r1-r0; % rho_w

cmap = color_scheme(length(k_all));
mn = 100;

% plot the simulated rho_w for a specific mn
figure;
hold on;
for iii = 1:length(k_all)
    id = find(m.*n == mn & k == k_all(iii));
    plot(log(m(id)./n(id))', r_w(id)','-o','Color',cmap(iii,:),'MarkerFaceColor',cmap(iii,:),...
        'LineWidth', 2, 'DisplayName', ['k = ', num2str(k_all(iii))]);
end
legend show
xlabel('log(m/n)');
ylabel('Simulated \rho_{diff}')
title(['mn = ', num2str(mn)])
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([-5 5 0 1])
box on;

% plot the fitted rho_w 
m = 2:0.1:mn/2;
n = mn./m;
r_min = (2*m+2*n-4)./(m.*n);
f0 = @(k) (r_min-1./(m.*n)) + (1-(r_min-1./(m.*n))).*exp(-6.3*sqrt(min(k./sqrt(m.*n),sqrt(m.*n))));
f1 = @(k) (r_min)  + (1-r_min).*exp(-0.15*min(k,m.*n).^(1.2)./sqrt(sqrt(m.*n)));
f = @(k) f1(k)-f0(k);
figure;
hold on;
for iii = 1:length(k_all)
    k = k_all(iii);
    plot(log(m./n)', f(k)','-','Color',cmap(iii,:),'LineWidth', 2, 'DisplayName', ['k = ', num2str(k)]);
end
legend show
xlabel('log(m/n)');
ylabel('Fitted \rho_w')
title(['mn = ', num2str(mn)])
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([-5 5 0 1])
box on;
