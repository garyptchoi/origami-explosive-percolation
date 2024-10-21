% Rigidity matrix building function with no planarity constraint
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

function A = RigidityMatrix(X, Y, Z, m, n)
    % Inputs:
    % X, Y, Z: Matrices containing the x, y, z coordinates of the grid vertices
    % m: Number of quads in the y-direction (vertical)
    % n: Number of quads in the x-direction (horizontal)
    
    M = m + 1; % Number of vertices in y-direction
    N = n + 1; % Number of vertices in x-direction

    %% Build the matrix G for horizontal edges
    
    row_all = [];
    col_all = [];
    val_all = [];
    
    for i = 1:n
        for j = 1:M  
            
            a = (j - 1) * n + i;  
    
            x_diff = 2 * (X(j, i) - X(j, i+1));
            y_diff = 2 * (Y(j, i) - Y(j, i+1));
            z_diff = 2 * (Z(j, i) - Z(j, i+1));
    
            col_indices = [
                (j - 1) * N * 3 + 3 * (i - 1) + 1,  
                (j - 1) * N * 3 + 3 * (i - 1) + 2,  
                (j - 1) * N * 3 + 3 * (i - 1) + 3,  
                (j - 1) * N * 3 + 3 * i + 1,        
                (j - 1) * N * 3 + 3 * i + 2,        
                (j - 1) * N * 3 + 3 * i + 3         
            ];
    
            
            values = [x_diff, y_diff, z_diff, -x_diff, -y_diff, -z_diff];
    
            row_all = [row_all, repmat(a, 1, 6)];
            col_all = [col_all, col_indices];
            val_all = [val_all, values];
        end
    end
    
    G = sparse(row_all, col_all, val_all, M * n, 3 * M * N);
    
    %% Build the matrix K for vertical edges
    
    row_all = [];
    col_all = [];
    val_all = [];   
    
    for i = 1:N  
        for j = 1:m  
    
            a = (j - 1) * N + i;  
    
            x_diff = 2 * (X(j, i) - X(j+1, i));
            y_diff = 2 * (Y(j, i) - Y(j+1, i));
            z_diff = 2 * (Z(j, i) - Z(j+1, i));
    
            col_indices = [
                (j - 1) * N * 3 + 3 * (i - 1) + 1,  
                (j - 1) * N * 3 + 3 * (i - 1) + 2,  
                (j - 1) * N * 3 + 3 * (i - 1) + 3,  
                j * N * 3 + 3 * (i - 1) + 1,        
                j * N * 3 + 3 * (i - 1) + 2,        
                j * N * 3 + 3 * (i - 1) + 3       
            ];
    
            
            values = [x_diff, y_diff, z_diff, -x_diff, -y_diff, -z_diff];
    
            row_all = [row_all, repmat(a, 1, 6)];
            col_all = [col_all, col_indices];
            val_all = [val_all, values];
        end
    end
    
    K = sparse(row_all, col_all, val_all, m * N, 3 * M * N);
   
    %% Build the matrix D for diagonal edges (no-shear constraints)
    
    row_all = [];
    col_all = [];
    val_all = [];
    
    for i = 1:n  
        for j = 1:m  
    
            a = (j - 1) * N + i;  
    
            x_diff = 2 * (X(j, i) - X(j+1, i+1));
            y_diff = 2 * (Y(j, i) - Y(j+1, i+1));
            z_diff = 2 * (Z(j, i) - Z(j+1, i+1));
    
            col_indices = [...
                (j - 1) * N * 3 + 3 * (i - 1) + 1,  
                (j - 1) * N * 3 + 3 * (i - 1) + 2,  
                (j - 1) * N * 3 + 3 * (i - 1) + 3,  
                j * N * 3 + 3 * i + 1,              
                j * N * 3 + 3 * i + 2,              
                j * N * 3 + 3 * i + 3               
            ];
    
            
            values = [x_diff, y_diff, z_diff, -x_diff, -y_diff, -z_diff];
    
            row_all = [row_all, repmat(a, 1, 6)];
            col_all = [col_all, col_indices];
            val_all = [val_all, values];
        end
    end
    
    D = sparse(row_all, col_all, val_all, M * N, 3 * M * N);
    
    %% Combine them
    A = [G; K; D];

end