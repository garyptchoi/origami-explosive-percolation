% Main program for running the simulation for the explosive rigidity 
% percolation in mxn Miura-ori structures (the rectangular case)
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

% pattern size (m*n)
mn_all = [25, 36, 49, 64, 81, 100, 121, 144, 169, 196, ...
    225, 256, 289, 324, 361, 400]; 

% number of choices in each step
k_all = [1, 2, 4, 8, 16, 32]; 

% selection rule (1: most efficient, 2: least efficient)
rule_all = [1,2]; 

% number of simulations for each setup
n_sim = 500;

for mn = mn_all
    mkdir(['simulation_results/rec_',num2str(mn)]);
end

%% Run the simulations (may take several hours for large mn and large k)

for rule = rule_all
    for mn = mn_all

        num_all = 2:mn-1;
        % get all combinations of m*n with m,n>=2
        L_all = num_all(rem(mn,num_all)==0); 

        for L = L_all
            %% Generate the Miura-ori structure
           
            m = L; % number of quads in y-direction
            n = mn/L; % number of quads in x-direction
            M = m+1; % number of vertices in y-direction
            N = n+1; % number of vertices in x-direction
            
            [X,Y,Z] = generate_miura_ori(M,N);
    
            num_quads = m*n;
           
            %% Construct the rigidity matrix of floppy miura origami
    
            A_initial = RigidityMatrix(X, Y, Z, m, n);
            DoF_initial = 3*M*N-calc_rank(A_initial)-6;
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Run simulations for each choice of k
            for k = k_all
               
                dof_all = zeros(n_sim,num_quads);
                
                tic;
                % using the parallel computing toolbox
                parfor t = 1:n_sim
                    disp(['Rule = ', num2str(rule), ', L = ', num2str(L),...
                        ', k = ', num2str(k),', Sim# = ',num2str(t)]);
                    
                    H = A_initial;
                    available_quads = true(num_quads, 1);
                
                    %% Iteratively add planarity constraints
                    
                    for q = 1:num_quads
                        % DoF of each selection at step Select_K
                        DoF_process=[];
                        
                        %% Step 1, Select_K 
                        available_idx = find(available_quads);
                    
                        % (New) if there are < K choices, still consider all of them and
                        % pick one based on the selection rule
                        quad_id_all = available_idx(randperm(length(available_idx), min(k,length(available_idx))));  
                            
                        j_all = floor((quad_id_all-1)/n)+1;
                        i_all = quad_id_all-(j_all-1)*n;
                
                        %% Step 2: Choose one among the k candidates
                        for p = 1:min(k,length(available_idx))
                        
                            i = i_all(p);
                            j = j_all(p);
                            
                            % update the rigidity matrix temporarily
                            V = computeV(i, j, X, Y, Z, M, N);                   
                            A = [H;V];
                            DoF = 3*M*N-calc_rank(A)-6;
                            DoF = min([DoF,DoF_initial]);
                            DoF = max([DoF,1]);
                            
                            DoF_process = [DoF_process,DoF];    
                        
                        end
                        
                        if rule == 1
                            % choose one that minimizes the DOF
                            [minValue, min_sel] = min(DoF_process);   
                            j= j_all(min_sel);
                            i = i_all(min_sel);
                        else
                            % choose one that maximizes the DOF
                            [maxValue, max_sel] = max(DoF_process);   
                            j= j_all(max_sel);
                            i = i_all(max_sel);
                        end
                    
                        % update the rigidity matrix based on the selection
                        V = computeV(i, j, X, Y, Z, M, N);
                        A = [H;V];
                        H = A;
                        DoF = 3*M*N-calc_rank(A)-6;
                        DoF = min([DoF,DoF_initial]);
                        DoF = max([DoF,1]);
                            
                        
                        dof_all(t,q) = DoF;
                        
                        selected_quad = (j-1)*n + i;  
                        available_quads(selected_quad) = false; 
                    
                    
                    end
            
                end
                
                % save result
                save_sim_rect(mn, L, k, rule, dof_all, n_sim);
            end
        end
    end
end
