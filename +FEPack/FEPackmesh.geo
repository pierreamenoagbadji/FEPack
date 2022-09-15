Include "/home/pierre/Documents/Etudes/these/codes/Chantier matlab/FEPack/+FEPack/+tools/FEPackGmsh_macros.geo";

h0 = 0.1;
is_structured = 0;
x1 = -2.00000000; y1 = 0.00000000; z1 = 0.0;
x2 = 0.50000000; y2 = 0.00000000; z2 = 0.0;
x3 = 0.50000000; y3 = 1.00000000; z3 = 0.0;
x4 = -2.00000000; y4 = 1.00000000; z4 = 0.0;
numNodesX = 32; numNodesY = 32;

domain_name = "rect";
side_name1 = "ymin";
side_name2 = "xmax";
side_name3 = "ymax";
side_name4 = "xmin";

Call FEPack_Rectangle;

Physical Point(1) = domain_0[];
Physical Line("ymin") = domain_1[];
Physical Line("xmax") = domain_2[];
Physical Line("ymax") = domain_3[];
Physical Line("xmin") = domain_4[];
Physical Surface("rect")= domain_5[];

Mesh.Format = 50;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;