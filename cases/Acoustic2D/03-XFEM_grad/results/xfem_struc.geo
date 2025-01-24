// Parameters: acoustic cavity
 
ax = 0.6;
hy = 0.6;

n1 = hy/20;

Point(1) = {ax,  0  , 0, n1};
Point(2) = {ax,  hy  , 0, n1};

Line(1) = {1, 2};

Physical Line(1) = {1};
Physical Point(2) = {2};
