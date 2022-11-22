Include "/home/pierre/Documents/Etudes/these/codes/Chantier matlab/FEPack/+FEPack/+tools/FEPackGmsh_macros.geo";

h0 = 0.1;
is_structured = 0;
x1 = -1.00000000; y1 = 0.00000000; z1 = 0.00000000;
x2 = 0.50000000; y2 = 0.00000000; z2 = 0.00000000;
x3 = 0.50000000; y3 = 1.00000000; z3 = 0.00000000;
x4 = -1.00000000; y4 = 1.00000000; z4 = 0.00000000;
x5 = -1.00000000; y5 = 0.00000000; z5 = 1.00000000;
x6 = 0.50000000; y6 = 0.00000000; z6 = 1.00000000;
x7 = 0.50000000; y7 = 1.00000000; z7 = 1.00000000;
x8 = -1.00000000; y8 = 1.00000000; z8 = 1.00000000;
numNodesX = 24; numNodesY = 16; numNodesZ = 16;

domain_name = "rect";
side_name1 = "xmax";
side_name2 = "xmin";
side_name3 = "ymax";
side_name4 = "ymin";
side_name5 = "zmax";
side_name6 = "zmin";

Call FEPack_Cuboid;

Physical Surface("xmax") = domain_1[];
Physical Surface("xmin") = domain_2[];
Physical Surface("ymax") = domain_3[];
Physical Surface("ymin") = domain_4[];
Physical Surface("zmax") = domain_5[];
Physical Surface("zmin") = domain_6[];
Physical Volume("rect")= domain_7[];

Mesh.Format = 50;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;