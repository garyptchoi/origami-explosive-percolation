% Adding a row vector to the rigidity matrix when a new planarity constraint is added
% V = the new planarity constraint added to the rigidity matrix as a row vector
% (i,j) = the quad which the new planarity constraint added to
%
% If you use this code in your work, please cite the following paper:
%    R. Li and G. P. T. Choi,
%    "Explosive rigidity percolation in origami."
%    Preprint, arXiv:2410.13945, 2024.
% 
% Copyright (c) 2024, Rongxuan Li and Gary P. T. Choi
% 
% https://github.com/garyptchoi/origami-explosive-percolation

function V = computeV(i, j, X, Y, Z, M, N)

    % v(j+1,i) --- v(j+1,i+1)
    %   (4)           (3)
    %    |             |
    %    |             |
    % v(j,i)   --- v(j,i+1)
    %   (1)           (2)

    % Want: 
    % (v3-v1) \dot (v2-v1)x(v4-v1) = 0
    % Determinant form:
    % |x3-x1  y3-y1  z3-z1|
    % |x2-x1  y2-y1  z2-z1| = 0
    % |x4-x1  y4-y1  z4-z1|

    size=3 * M * N;
       
    y21 = Y(j, i+1) - Y(j, i);  
    y31 = Y(j+1, i+1) - Y(j, i);
    y41 = Y(j+1, i) - Y(j, i);
    
    z21 = Z(j, i+1) - Z(j, i);
    z31 = Z(j+1, i+1) - Z(j, i);
    z41 = Z(j+1, i) - Z(j, i);
    
    x21 = X(j, i+1) - X(j, i);  
    x31 = X(j+1, i+1) - X(j, i);
    x41 = X(j+1, i) - X(j, i);

    a=(j-1) * N * 3 + 3 * (i-1) + 1;
    b=(j-1) * N * 3 + 3 * (i) + 1;
    c=(j) * N * 3 + 3 * (i) + 1;
    d=(j) * N * 3 + 3 * (i-1) + 1;
    
    a1 = -(y21 * z41 - y41 * z21) + (y31 * z41 - y41 * z31) - (y31 * z21 - y21 * z31);
    a2 =  (x21 * z41 - x41 * z21) - (x31 * z41 - x41 * z31) + (x31 * z21 - x21 * z31);
    a3 = -(x21 * y41 - x41 * y21) + (x31 * y41 - x41 * y31) - (x31 * y21 - x21 * y31);
    
    b1 = -(y31 * z41 - y41 * z31);
    b2 =  x31 * z41 - x41 * z31;
    b3 = -(x31 * y41 - x41 * y31);
    
    c1 =  y21 * z41 - y41 * z21;
    c2 = -(x21 * z41 - x41 * z21);
    c3 =  x21 * y41 - x41 * y21;
    
    d1 =  y31 * z21 - y21 * z31;
    d2 = -(x31 * z21 - x21 * z31);
    d3 =  x31 * y21 - x21 * y31;

    row=[1,1,1,1,1,1,1,1,1,1,1,1];
    col=[a,a+1,a+2,b,b+1,b+2,c,c+1,c+2,d,d+1,d+2];
    val=[a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3];

    V=sparse(row, col, val, 1, size);


end
