% Analysis of the simulation results with different geometric parameters
% 
% Procedure:
% 1. Changing the parameters gamma and theta in generate_miura_ori.m
%    Alternative setup 1: (gamma,theta) = (pi/3, pi/3) 
%    (results saved in square_pi3)
%    Alternative setup 2: (gamma,theta) = (pi/6, pi/6)
%    (results saved in square_pi6)
%
% 2. Load the simulation results and analyze them
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
L_all = 10;

% number of choices in each step
k_all = [1, 2, 4, 8, 16, 32];

% selection rule (1: most efficient, 2: least efficient)
rule_all = [1, 2]; 

% number of simulations for each setup
n_sim = 500;

%% Plot the Miura-ori structure

% gamma = pi/4; theta = acos(sqrt(2/3)); % default
gamma = pi/3; theta = pi/3; % alternative setup 1
% gamma = pi/6; theta = pi/6; % alternative setup 2

L = 10;
m = L;
n = L;
M = m+1;
N = n+1;
[X,Y,Z] = generate_miura_ori(M,N,gamma,theta);

figure;
hold on;
for i = 1:m
    for j = 1:n

        x_quad = [X(i,j), X(i,j+1), X(i+1,j+1), X(i+1,j)];
        y_quad = [Y(i,j), Y(i,j+1), Y(i+1,j+1), Y(i+1,j)];
        z_quad = [Z(i,j), Z(i,j+1), Z(i+1,j+1), Z(i+1,j)];

        fill3(x_quad, y_quad, z_quad, [255 230 213]/255, ...
            'EdgeColor','k','LineWidth',2);

    end
end

axis equal off; 
ax = gca; ax.Clipping = 'off';
view([24 22])

%% Plot r vs P for each fixed rule and each fixed L

cmap = color_scheme(length(k_all));

for rule = rule_all
    for L = L_all
        figure;
        hold on;
        for k = k_all
        
            id = find(k_all == k);
            
            % choose one set of simulation results to load
%             load(['simulation_results/square/sim_L_',num2str(L),'_k_',num2str(k),'_rule_', num2str(rule),'.mat']);
            load(['simulation_results/square_pi3/sim_L_',num2str(L),'_k_',num2str(k),'_rule_', num2str(rule),'.mat']);
%             load(['simulation_results/square_pi6/sim_L_',num2str(L),'_k_',num2str(k),'_rule_', num2str(rule),'.mat']);

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

    end
end




