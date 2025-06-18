# Explosive rigidity percolation in origami

<img src = "https://github.com/garyptchoi/origami-explosive-percolation/blob/main/cover.jpg" height="300" />

This repository contains simulation codes for rigidifying a floppy Miura-ori structure by incrementally adding facet planarity constraints.
At every step, k random candidate facets are chosen and one among them is selected based on some selection rule.
Changing the number of choices k and the selection rule allows us to effectively control the explosive rigidity percolation behavior in origami.

Any comments and suggestions are welcome. 

If you use this code in your work, please cite the following paper:

R. Li and G. P. T. Choi,
"[Explosive rigidity percolation in origami.](https://doi.org/10.1098/rspa.2024.0826)"
Proceedings of the Royal Society A, 481(2316), 20240826, 2025. 

Copyright (c) 2024-2025, Rongxuan Li and Gary P. T. Choi

https://github.com/garyptchoi/origami-explosive-percolation

===============================================================

Simulation codes:
* `run_simulation_sqaure.m`: The main program for the simulation for the square case ($L \times L$ Miura-ori structure)
* `run_simulation_rect.m`: The main program for the simulation for the rectangular case ($m \times n$ Miura-ori structure)

Analysis codes:
* `analysis_sqaure.m`: The main program for the analysis of the square case ($L \times L$ Miura-ori structure)
* `analysis_rect.m`: The main program for the analysis of the rectangular case ($m \times n$ Miura-ori structure)
* `analysis_geom.m`: The main program for the analysis of alternative geometric parameters 

Simulation results:
* `simulation_results/square/sim_L_aaa_k_bbb_rule_ccc.mat`: The simulation result for the square case (with pattern size L = aaa, number of choices k = bbb, rule = ccc)
* `simulation_results/rec_ddd/sim_rect_ddd_L_eee_k_fff_rule_ggg.mat': The simulation result for the rectangular case (with pattern size S = ddd = $m \times n$ where m = eee, number of choices k = fff, rule = ggg)
* `simulation_results/square_pi3`, `simulation_results/square_pi6`: the simulation results for Miura-ori structures with alternative geometric parameters
