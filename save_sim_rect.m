% Save the rectangular simulation results
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

function save_sim_rect(mn,L,k,rule,dof_all,n_sim) 
    save(['simulation_results/rec_',num2str(mn),'/sim_rect_',num2str(mn),...
        '_L_',num2str(L),'_k_',num2str(k),'_rule_',num2str(rule),'.mat'],...
        'mn','L','dof_all','k','n_sim');
end
