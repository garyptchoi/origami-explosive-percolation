% Save the simulation results
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

function save_sim(L,k,rule,dof_all,n_sim) 
    save(['simulation_results/square/sim_L_',num2str(L),'_k_',num2str(k),...
        '_rule_',num2str(rule),'.mat'],'L','dof_all','k','n_sim');
end
