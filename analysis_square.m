% Analysis of the simulation results for the explosive rigidity percolation
% in LxL Miura-ori structures (the square case)
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

% pattern size
L_all = [5, 10, 15, 20, 25, 30];

% number of choices in each step
k_all = [1, 2, 4, 8, 16, 32];

% selection rule (1: most efficient, 2: least efficient)
rule_all = [1,2]; 

% number of simulations for each setup
n_sim = 500;

%% Plot r vs (d-1)/(d_initial-1) for each rule and each k

cmap = color_scheme(length(k_all));

for rule = rule_all
    for k = k_all
        figure;
        hold on;
        for L = L_all
            
            id = find(L_all == L);
    
            load(['simulation_results/square/sim_L_',num2str(L),'_k_',...
                num2str(k),'_rule_', num2str(rule),'.mat']);

            m = L;
            n = L;
            M = m+1;
            N = n+1;
            [X,Y,Z] = generate_miura_ori(M,N);
            A_initial = RigidityMatrix(X, Y, Z, m, n);
            DoF_initial = 3*M*N-calc_rank(A_initial)-6;
            
            plot(repmat(linspace(0,1,m*n+1),n_sim,1)',...
                [DoF_initial*ones(n_sim,1)-1,dof_all-1]'/(DoF_initial-1),...
                'Color',[cmap(id,:),0.1], 'LineWidth',1);
        end
        title(['k = ', num2str(k)]);
        xlabel('\rho');
        ylabel('(d-1)/(d_{initial}-1)');
        set(gca,'FontSize',20);
        set(gca,'LineWidth',2);
    end
end

%% Plot r vs P for each fixed rule and each fixed L

cmap = color_scheme(length(k_all));

for rule = rule_all
    for L = L_all
        figure;
        hold on;
        for k = k_all
        
            id = find(k_all == k);
        
            load(['simulation_results/square/sim_L_',num2str(L),'_k_',...
                num2str(k),'_rule_', num2str(rule),'.mat']);

            m = L;
            n = L;
            M = m+1;
            N = n+1;
            [X,Y,Z] = generate_miura_ori(M,N);
            A_initial = RigidityMatrix(X, Y, Z, m, n);
            DoF_initial = 3*M*N-calc_rank(A_initial)-6;
            
            % planarity constraint density
            r = linspace(0,1,m*n+1)';
            
            % probability of getting a 1-DoF structure
            P = sum([DoF_initial*ones(n_sim,1),dof_all]==1,1)/n_sim;
        
            plot(r,P,'Color',cmap(id,:),'LineWidth',3,'DisplayName',['k = ',num2str(k)]);
        
        end
       
        legend show;
        title(['L = ', num2str(L)]);
        xlabel('\rho');
        ylabel('P');
        set(gca,'FontSize',20);
        set(gca,'LineWidth',2);
        axis([0 1 0 1])
    end
end


%% Plot k vs rho^* for each fixed rule and each fixed L

LL = repmat(L_all', 1, length(k_all));
kk = repmat(k_all, length(L_all), 1);
kk(kk>=LL.^2) = (LL(kk>=LL.^2)).^2;

r_star_all = zeros(length(L_all),length(k_all)); % rho^*

r_0_all = zeros(length(L_all),length(k_all)); % rho^0
r_1_all = zeros(length(L_all),length(k_all)); % rho_1

r_01_all = zeros(length(L_all),length(k_all)); % rho^0.1
r_09_all = zeros(length(L_all),length(k_all)); % rho_0.9

r_025_all = zeros(length(L_all),length(k_all)); % rho^0.25
r_075_all = zeros(length(L_all),length(k_all)); % rho_0.75

for rule = [2,1] % run rule = 2 first for the plots, as the rule = 1 results will be further analyzed below
    for L = L_all
        id_L = find(L_all == L);
        for k = k_all
        
            id_k = find(k_all == k);
        
            load(['simulation_results/square/sim_L_',num2str(L),'_k_',...
                num2str(k),'_rule_', num2str(rule),'.mat']);

            m = L;
            n = L;
            M = m+1;
            N = n+1;
            [X,Y,Z] = generate_miura_ori(M,N);
            A_initial = RigidityMatrix(X, Y, Z, m, n);
            DoF_initial = 3*M*N - calc_rank(A_initial) - 6;

            % percentage of planarity constraints added
            r = linspace(0,1,m*n+1)';
            
            % probability of getting a 1-DoF structure
            P = sum([DoF_initial*ones(n_sim,1),dof_all]==1,1)/n_sim;
            
            % critical rho
            idx = find(P >= 0.5, 1); 
            r_star_all(id_L,id_k) = r(idx);

            % r^0, r_1
            idx_r0 = find(P == 0,1,'last'); 
            idx_r1 = find(P == 1, 1); 
            r_0_all(id_L,id_k) = r(idx_r0);
            r_1_all(id_L,id_k) = r(idx_r1);

            % r^0.1, r_0.9
            id_r1 = find(P <= 0.1, 1,'last');  
            id_r9 = find(P >= 0.9, 1); 
            r_01_all(id_L,id_k) = r(id_r1);
            r_09_all(id_L,id_k) = r(id_r9);

            % r^0.25, r_0.75
            id_r25 = find(P <= 0.25, 1,'last');   
            id_r75 = find(P >= 0.75, 1); 
            r_025_all(id_L,id_k) = r(id_r25);
            r_075_all(id_L,id_k) = r(id_r75);

        end

        % Plot k vs critical rho
        figure;
        if rule == 1
            plot([0,34],[(4*L-4)/L^2,(4*L-4)/L^2],'k--','LineWidth',1);
            hold on;
            xlim([0 34])
            ylim([0 1])
            plot(k_all, r_star_all(id_L,:), '-o', ...
                'Color', [201 0 22]/255, 'MarkerFaceColor', [201 0 22]/255, ...
                'LineWidth', 2);
            xlabel('k');
            ylabel('\rho*')
            title(['L = ', num2str(L)]);
            set(gca,'FontSize',18);
            set(gca,'LineWidth',2);

        end

        if rule == 2
            plot(k_all, r_star_all(id_L,:), '-o', ...
                'Color', [201 0 22]/255, 'MarkerFaceColor', [201 0 22]/255, ...
                'LineWidth', 2);
            xlim([0 34])
            ylim([0 1])
            xlabel('k');
            ylabel('\rho*')
            title(['L = ', num2str(L)]);
            set(gca,'FontSize',18);
            set(gca,'LineWidth',2);
        end

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyzing rho_diff

% plot the simulated rho_diff vs k for each L

cmap = color_scheme(length(L_all));

r_diff = r_star_all - ((4*L_all-4)./L_all.^2)';

figure;
hold on;
for i = 1:size(r_diff,1)
    L = L_all(i);
    plot(k_all', r_diff(i,:),'-o','Color',cmap(i,:),...
        'MarkerFaceColor',cmap(i,:),'LineWidth', 2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Simulated \rho_{diff}')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

% plot the fitted rho_diff vs k for each L

k_range = 1:0.1:32;

figure;
hold on;
for i = 1:size(r_diff,1)
    L = L_all(i);
    rho_diff_fit = (1+min(k_range/L,L)*7.8).^(-1.4);
    plot(k_range', rho_diff_fit','-','Color',cmap(i,:),...
        'LineWidth', 2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Fitted \rho_{diff}')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyzing rho_w

% plot the simulated rho_w vs k for each L

cmap = color_scheme(length(L_all));

r_w = r_1_all - r_0_all;

figure;
hold on;
for i = 1:size(r_0_all,1)
    L = L_all(i);
    plot(k_all', r_0_all(i,:),'-o','Color',cmap(i,:),...
        'MarkerFaceColor',cmap(i,:),'LineWidth', 2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Simulated \rho^{0}')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

figure;
hold on;
for i = 1:size(r_1_all,1)
    L = L_all(i);
    plot(k_all', r_1_all(i,:),'-o','Color',cmap(i,:),...
        'MarkerFaceColor',cmap(i,:),'LineWidth', 2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Simulated \rho_{1}')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

figure;
hold on;
for i = 1:size(r_w,1)
    L = L_all(i);
    plot((k_all)', r_w(i,:),'-o','Color',cmap(i,:),...
        'MarkerFaceColor',cmap(i,:),'LineWidth', 2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Simulated \rho_{w}')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

% plot the fitted model
k_range = 1:0.1:32;

figure;
hold on;
for iii = 1:length(L_all)
    L = L_all(iii);
    rmin = (4*L-4)/L^2;
    fit_r0 = (rmin-1/L^2)  + (1-(rmin-1/L^2))*exp(-6.3*min(k_range,L^2).^(0.5)./L.^(0.5));

    plot(k_range', fit_r0','-','Color',cmap(iii,:),'LineWidth',2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Fitted \rho^0')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

figure;
hold on;
for iii = 1:length(L_all)
    L = L_all(iii);
    rmin = (4*L-4)/L^2;
    fit_r1 = (rmin)  + (1-rmin)*exp(-0.15*min(k_range,L^2).^(1.2)./sqrt(L));

    plot(k_range', fit_r1','-','Color',cmap(iii,:),'LineWidth',2,'DisplayName',['L = ',num2str(L)]);
end
legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Fitted \rho_1')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;

figure;
hold on;
for iii = 1:length(L_all)
    L = L_all(iii);
    rmin = (4*L-4)/L^2;
    fit_r0 = (rmin-1/L^2) + (1-(rmin-1/L^2))*exp(-6.3*sqrt(min(k_range/L,L)));
    fit_r1 = (rmin) + (1-rmin)*exp(-0.15*min(k_range,L^2).^(1.2)./sqrt(L));
    fit_rw = fit_r1-fit_r0;

    plot(k_range', fit_rw','-','Color',cmap(iii,:),'LineWidth',2,'DisplayName',['L = ',num2str(L)]);
end

legend show;
xlim([0 35])
ylim([0 1])
xlabel('k');
ylabel('Fitted \rho_w')
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
box on;
