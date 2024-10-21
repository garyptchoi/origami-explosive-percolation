% Miura-ori structure generator
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

function [X,Y,Z] = generate_miura_ori(M,N,gamma,theta)
        a = 2; % length of parallelogram side
        b = 2; % length of the other side

        if nargin < 3
            gamma = pi/4; % angle between sides
        end
        if nargin < 4
            theta = acos(sqrt(2/3)); % folding angle
        end
        
        X = zeros(M, N);
        Y = zeros(M, N);
        Z = zeros(M, N);
        
        for i = 1:M
            for j = 1:N
                X(i,j) = (j-1) * a * cos(gamma);
             if mod(j, 2) == 0    
                Y(i,j) = (i-1) * b; 
             else
                Y(i,j) = (i-1) * b +  a * sin(gamma);
             end
                
                if mod(i, 2) == 0
                    Z(i,j) = 0; 
                else
                    Z(i,j) = a * sin(theta); 
                end
            end
        end
end