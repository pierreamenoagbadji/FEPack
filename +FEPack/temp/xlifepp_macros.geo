General.ExpertMode=1;

// global variables
itp=0; // iterator on points
itc=0; // iterator on curves
its=0; // iterator on surfaces
nbp=1000; // maximum number of points
nbe=2000; // maximum number of elliptic arcs
nbl=2000; // maximum number of lines
nbs=1000; // maximum number of surfaces
nbloops=6000; // maximum number of loops
nbout=2000; // maximum size of array returned by extrusion
h0=0.1; // default characterictic length

// every builtX_i variable is initialized to -1
Function xlifepp_init
  For i In {1:nbp}
    builtP~{i}=-1;
  EndFor
  For i In {1:nbe}
    builtE~{i}=-1;
  EndFor
  For i In {1:nbl}
    builtL~{i}=-1;
    builtC~{i}=-1;
  EndFor
  For i In {1:nbs}
    builtS~{i}=-1;
  EndFor
  For i In {0:nbloops}
    builtloops~{i}=-1;
  EndFor
  For i In {0:nbout}
    builtout[i]=-1;
    builtouta[i]=-1;
    builtoutb[i]=-1;
    builtoutc[i]=-1;
  EndFor
  LoopsStorage=1;
  Apogee=-1;
Return

Function xlifepp_Segment
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If(builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  If (L_1 < 0)
    L_1=-L_1;
  EndIf
  If (LoopsStorage == 1)
    loops~{l}=L_1;
  EndIf
  curves~{itc}=L_1;
  itc++;
Return

Function xlifepp_EllArc
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  If (builtP_3 == -1)
    If (Apogee == 1)
      P_3=P_1;
    EndIf
    If (Apogee != 1)
      P_3=newp;
      Point(P_3)={x3,y3,z3,h3};
    EndIf
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  If (builtP_4 == -1)
    If (Apogee == 2)
      P_4=P_3;
    EndIf
    If (Apogee != 2)
      P_4=newp;
      Point(P_4)={x4,y4,z4,h4};
    EndIf
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_1,P_2,P_3,P_4};
  EndIf
  If(builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  If (E_1 < 0)
    E_1=-E_1;
  EndIf
  If (LoopsStorage == 1)
    loops~{l}=E_1;
  EndIf
  curves~{itc}=E_1;
  itc++;
Return

Function xlifepp_CircArc
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtE_1 == -1)
    E_1=newl;
    Circle(E_1)={P_1,P_2,P_3};
  EndIf
  If(builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  If (E_1 < 0)
    E_1=-E_1;
  EndIf
  If (LoopsStorage == 1)
    loops~{l}=E_1;
  EndIf
  curves~{itc}=E_1;
  itc++;
Return

Function xlifepp_Triangle
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_2,P_3};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_3,P_1};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;

  LL_1=newll;
  If (LoopsStorage ==1)
    loops~{l}=LL_1;
  EndIf
  Line Loop(LL_1)={L_1,L_2,L_3};
  surfs~{its}=LL_1;
  its++;
Return

Function xlifepp_Quadrangle
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_2,P_3};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_3,P_4};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_4,P_1};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  If (LoopsStorage == 1)
    loops~{l}=LL_1;
  EndIf
  Line Loop(LL_1)={L_1,L_2,L_3,L_4};
  surfs~{its}=LL_1;
  its++;
Return

Function xlifepp_Ellipse
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;

  LL_1=newll;
  If (LoopsStorage == 1)
    loops~{l}=LL_1;
  EndIf
  Line Loop(LL_1)={E_1,E_2,E_3,E_4};
  surfs~{its}=LL_1;
  its++;
Return

Function xlifepp_AcuteEllipticSector
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_3};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_2,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_4,P_1};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;

  LL_1=newll;
  If (LoopsStorage == 1)
    loops~{l}=LL_1;
  EndIf
  Line Loop(LL_1)={L_1,E_2,L_3};
  surfs~{its}=LL_1;
  its++;
Return

Function xlifepp_ObtuseEllipticSector
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_3};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_2,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_2,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_5,P_1};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  If (LoopsStorage == 1)
    loops~{l}=LL_1;
  EndIf
  Line Loop(LL_1)={L_1,E_2,E_3,L_4};
  surfs~{its}=LL_1;
  its++;
Return

Function xlifepp_Tetrahedron
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_2,P_3};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_3,P_1};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_1,P_4};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_2,P_4};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;
  If (builtL_6 == -1)
    L_6=newl;
    Line(L_6)={P_3,P_4};
  EndIf
  If (builtL_6 != -1)
    L_6=builtL_6;
  EndIf
  curves~{itc}=L_6;
  If (curves~{itc} < 0)
    curves~{itc}=-L_6;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-L_1,-L_2,-L_3};
  LL_2=newll;
  Line Loop(LL_2)={L_1,L_5,-L_4};
  LL_3=newll;
  Line Loop(LL_3)={L_2,L_6,-L_5};
  LL_4=newll;
  Line Loop(LL_4)={L_3,L_4,-L_6};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4};
Return

Function xlifepp_Hexahedron
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_2,P_3};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_3,P_4};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_4,P_1};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_5,P_6};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;
  If (builtL_6 == -1)
    L_6=newl;
    Line(L_6)={P_6,P_7};
  EndIf
  If (builtL_6 != -1)
    L_6=builtL_6;
  EndIf
  curves~{itc}=L_6;
  If (curves~{itc} < 0)
    curves~{itc}=-L_6;
  EndIf
  itc++;
  If (builtL_7 == -1)
    L_7=newl;
    Line(L_7)={P_7,P_8};
  EndIf
  If (builtL_7 != -1)
    L_7=builtL_7;
  EndIf
  curves~{itc}=L_7;
  If (curves~{itc} < 0)
    curves~{itc}=-L_7;
  EndIf
  itc++;
  If (builtL_8 == -1)
    L_8=newl;
    Line(L_8)={P_8,P_5};
  EndIf
  If (builtL_8 != -1)
    L_8=builtL_8;
  EndIf
  curves~{itc}=L_8;
  If (curves~{itc} < 0)
    curves~{itc}=-L_8;
  EndIf
  itc++;
  If (builtL_9 == -1)
    L_9=newl;
    Line(L_9)={P_1,P_5};
  EndIf
  If (builtL_9 != -1)
    L_9=builtL_9;
  EndIf
  curves~{itc}=L_9;
  If (curves~{itc} < 0)
    curves~{itc}=-L_9;
  EndIf
  itc++;
  If (builtL_10 == -1)
    L_10=newl;
    Line(L_10)={P_2,P_6};
  EndIf
  If (builtL_10 != -1)
    L_10=builtL_10;
  EndIf
  curves~{itc}=L_10;
  If (curves~{itc} < 0)
    curves~{itc}=-L_10;
  EndIf
  itc++;
  If (builtL_11 == -1)
    L_11=newl;
    Line(L_11)={P_3,P_7};
  EndIf
  If (builtL_11 != -1)
    L_11=builtL_11;
  EndIf
  curves~{itc}=L_11;
  If (curves~{itc} < 0)
    curves~{itc}=-L_11;
  EndIf
  itc++;
  If (builtL_12 == -1)
    L_12=newl;
    Line(L_12)={P_4,P_8};
  EndIf
  If (builtL_12 != -1)
    L_12=builtL_12;
  EndIf
  curves~{itc}=L_12;
  If (curves~{itc} < 0)
    curves~{itc}=-L_12;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-L_1,-L_2,-L_3,-L_4};
  LL_2=newll;
  Line Loop(LL_2)={L_5,L_6,L_7,L_8};
  LL_3=newll;
  Line Loop(LL_3)={L_1,L_10,-L_5,-L_9};
  LL_4=newll;
  Line Loop(LL_4)={L_3,L_12,-L_7,-L_11};
  LL_5=newll;
  Line Loop(LL_5)={L_4,L_9,-L_8,-L_12};
  LL_6=newll;
  Line Loop(LL_6)={L_2,L_11,-L_6,-L_10};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6};
Return

Function xlifepp_Parallelepiped7
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;
  If (builtP_9 == -1)
    P_9=newp;
    Point(P_9)={x9,y9,z9,h9};
  EndIf
  If (builtP_9 != -1)
    P_9=builtP_9;
  EndIf
  points~{itp}=P_9;
  itp++;
  If (builtP_10 == -1)
    P_10=newp;
    Point(P_10)={x10,y10,z10,h10};
  EndIf
  If (builtP_10 != -1)
    P_10=builtP_10;
  EndIf
  points~{itp}=P_10;
  itp++;
  If (builtP_11 == -1)
    P_11=newp;
    Point(P_11)={x11,y11,z11,h11};
  EndIf
  If (builtP_11 != -1)
    P_11=builtP_11;
  EndIf
  points~{itp}=P_11;
  itp++;
  If (builtP_12 == -1)
    P_12=newp;
    Point(P_12)={x12,y12,z12,h12};
  EndIf
  If (builtP_12 != -1)
    P_12=builtP_12;
  EndIf
  points~{itp}=P_12;
  itp++;
  If (builtP_13 == -1)
    P_13=newp;
    Point(P_13)={x13,y13,z13,h13};
  EndIf
  If (builtP_13 != -1)
    P_13=builtP_13;
  EndIf
  points~{itp}=P_13;
  itp++;
  If (builtP_14 == -1)
    P_14=newp;
    Point(P_14)={x14,y14,z14,h14};
  EndIf
  If (builtP_14 != -1)
    P_14=builtP_14;
  EndIf
  points~{itp}=P_14;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_8};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_8,P_11};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_11,P_9};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_9,P_2};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_2,P_3};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;
  If (builtL_6 == -1)
    L_6=newl;
    Line(L_6)={P_3,P_1};
  EndIf
  If (builtL_6 != -1)
    L_6=builtL_6;
  EndIf
  curves~{itc}=L_6;
  If (curves~{itc} < 0)
    curves~{itc}=-L_6;
  EndIf
  itc++;
  If (builtL_7 == -1)
    L_7=newl;
    Line(L_7)={P_4,P_5};
  EndIf
  If (builtL_7 != -1)
    L_7=builtL_7;
  EndIf
  curves~{itc}=L_7;
  If (curves~{itc} < 0)
    curves~{itc}=-L_7;
  EndIf
  itc++;
  If (builtL_8 == -1)
    L_8=newl;
    Line(L_8)={P_5,P_6};
  EndIf
  If (builtL_8 != -1)
    L_8=builtL_8;
  EndIf
  curves~{itc}=L_8;
  If (curves~{itc} < 0)
    curves~{itc}=-L_8;
  EndIf
  itc++;
  If (builtL_9 == -1)
    L_9=newl;
    Line(L_9)={P_6,P_7};
  EndIf
  If (builtL_9 != -1)
    L_9=builtL_9;
  EndIf
  curves~{itc}=L_9;
  If (curves~{itc} < 0)
    curves~{itc}=-L_9;
  EndIf
  itc++;
  If (builtL_10 == -1)
    L_10=newl;
    Line(L_10)={P_7,P_4};
  EndIf
  If (builtL_10 != -1)
    L_10=builtL_10;
  EndIf
  curves~{itc}=L_10;
  If (curves~{itc} < 0)
    curves~{itc}=-L_10;
  EndIf
  itc++;
  If (builtL_11 == -1)
    L_11=newl;
    Line(L_11)={P_12,P_10};
  EndIf
  If (builtL_11 != -1)
    L_11=builtL_11;
  EndIf
  curves~{itc}=L_11;
  If (curves~{itc} < 0)
    curves~{itc}=-L_11;
  EndIf
  itc++;
  If (builtL_12 == -1)
    L_12=newl;
    Line(L_12)={P_10,P_13};
  EndIf
  If (builtL_12 != -1)
    L_12=builtL_12;
  EndIf
  curves~{itc}=L_12;
  If (curves~{itc} < 0)
    curves~{itc}=-L_12;
  EndIf
  itc++;
  If (builtL_13 == -1)
    L_13=newl;
    Line(L_13)={P_13,P_14};
  EndIf
  If (builtL_13 != -1)
    L_13=builtL_13;
  EndIf
  curves~{itc}=L_13;
  If (curves~{itc} < 0)
    curves~{itc}=-L_13;
  EndIf
  itc++;
  If (builtL_14 == -1)
    L_14=newl;
    Line(L_14)={P_14,P_12};
  EndIf
  If (builtL_14 != -1)
    L_14=builtL_14;
  EndIf
  curves~{itc}=L_14;
  If (curves~{itc} < 0)
    curves~{itc}=-L_14;
  EndIf
  itc++;
  If (builtL_15 == -1)
    L_15=newl;
    Line(L_15)={P_1,P_4};
  EndIf
  If (builtL_15 != -1)
    L_15=builtL_15;
  EndIf
  curves~{itc}=L_15;
  If (curves~{itc} < 0)
    curves~{itc}=-L_15;
  EndIf
  itc++;
  If (builtL_16 == -1)
    L_16=newl;
    Line(L_16)={P_10,P_5};
  EndIf
  If (builtL_16 != -1)
    L_16=builtL_16;
  EndIf
  curves~{itc}=L_16;
  If (curves~{itc} < 0)
    curves~{itc}=-L_16;
  EndIf
  itc++;
  If (builtL_17 == -1)
    L_17=newl;
    Line(L_17)={P_2,P_6};
  EndIf
  If (builtL_17 != -1)
    L_17=builtL_17;
  EndIf
  curves~{itc}=L_17;
  If (curves~{itc} < 0)
    curves~{itc}=-L_17;
  EndIf
  itc++;
  If (builtL_18 == -1)
    L_18=newl;
    Line(L_18)={P_3,P_7};
  EndIf
  If (builtL_18 != -1)
    L_18=builtL_18;
  EndIf
  curves~{itc}=L_18;
  If (curves~{itc} < 0)
    curves~{itc}=-L_18;
  EndIf
  itc++;
  If (builtL_19 == -1)
    L_19=newl;
    Line(L_19)={P_8,P_12};
  EndIf
  If (builtL_19 != -1)
    L_19=builtL_19;
  EndIf
  curves~{itc}=L_19;
  If (curves~{itc} < 0)
    curves~{itc}=-L_19;
  EndIf
  itc++;
  If (builtL_20 == -1)
    L_20=newl;
    Line(L_20)={P_11,P_14};
  EndIf
  If (builtL_20 != -1)
    L_20=builtL_20;
  EndIf
  curves~{itc}=L_20;
  If (curves~{itc} < 0)
    curves~{itc}=-L_20;
  EndIf
  itc++;
  If (builtL_21 == -1)
    L_21=newl;
    Line(L_21)={P_9,P_13};
  EndIf
  If (builtL_21 != -1)
    L_21=builtL_21;
  EndIf
  curves~{itc}=L_21;
  If (curves~{itc} < 0)
    curves~{itc}=-L_21;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-L_1,-L_2,-L_3,-L_4,-L_5,-L_6};
  LL_2=newll;
  Line Loop(LL_2)={L_7,L_8,L_9,L_10};
  LL_3=newll;
  Line Loop(LL_3)={L_1,L_19,L_11,L_16,-L_7,-L_15};
  LL_4=newll;
  Line Loop(LL_4)={L_5,L_18,-L_9,-L_17};
  LL_5=newll;
  Line Loop(LL_5)={L_6,L_15,-L_10,-L_18};
  LL_6=newll;
  Line Loop(LL_6)={L_4,L_17,-L_8,-L_16,L_12,-L_21};
  LL_7=newll;
  Line Loop(LL_7)={-L_11,-L_12,-L_13,-L_14};
  LL_8=newll;
  Line Loop(LL_8)={L_3,L_21,L_13,-L_20};
  LL_9=newll;
  Line Loop(LL_9)={L_2,L_20,L_14,-L_19};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;
  If (builtS_9 == -1)
    S_9=news;
    Plane Surface(S_9)={LL_9};
  EndIf
  If (builtS_9 != -1)
    S_9=builtS_9;
  EndIf
  surfs~{its}=S_9;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8,S_9};
Return

Function xlifepp_Parallelepiped6
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;
  If (builtP_9 == -1)
    P_9=newp;
    Point(P_9)={x9,y9,z9,h9};
  EndIf
  If (builtP_9 != -1)
    P_9=builtP_9;
  EndIf
  points~{itp}=P_9;
  itp++;
  If (builtP_10 == -1)
    P_10=newp;
    Point(P_10)={x10,y10,z10,h10};
  EndIf
  If (builtP_10 != -1)
    P_10=builtP_10;
  EndIf
  points~{itp}=P_10;
  itp++;
  If (builtP_11 == -1)
    P_11=newp;
    Point(P_11)={x11,y11,z11,h11};
  EndIf
  If (builtP_11 != -1)
    P_11=builtP_11;
  EndIf
  points~{itp}=P_11;
  itp++;
  If (builtP_12 == -1)
    P_12=newp;
    Point(P_12)={x12,y12,z12,h12};
  EndIf
  If (builtP_12 != -1)
    P_12=builtP_12;
  EndIf
  points~{itp}=P_12;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_7,P_9};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_9,P_1};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_2};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_2,P_7};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_3,P_4};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;
  If (builtL_6 == -1)
    L_6=newl;
    Line(L_6)={P_4,P_5};
  EndIf
  If (builtL_6 != -1)
    L_6=builtL_6;
  EndIf
  curves~{itc}=L_6;
  If (curves~{itc} < 0)
    curves~{itc}=-L_6;
  EndIf
  itc++;
  If (builtL_7 == -1)
    L_7=newl;
    Line(L_7)={P_5,P_6};
  EndIf
  If (builtL_7 != -1)
    L_7=builtL_7;
  EndIf
  curves~{itc}=L_7;
  If (curves~{itc} < 0)
    curves~{itc}=-L_7;
  EndIf
  itc++;
  If (builtL_8 == -1)
    L_8=newl;
    Line(L_8)={P_6,P_3};
  EndIf
  If (builtL_8 != -1)
    L_8=builtL_8;
  EndIf
  curves~{itc}=L_8;
  If (curves~{itc} < 0)
    curves~{itc}=-L_8;
  EndIf
  itc++;
  If (builtL_9 == -1)
    L_9=newl;
    Line(L_9)={P_8,P_10};
  EndIf
  If (builtL_9 != -1)
    L_9=builtL_9;
  EndIf
  curves~{itc}=L_9;
  If (curves~{itc} < 0)
    curves~{itc}=-L_9;
  EndIf
  itc++;
  If (builtL_10 == -1)
    L_10=newl;
    Line(L_10)={P_10,P_12};
  EndIf
  If (builtL_10 != -1)
    L_10=builtL_10;
  EndIf
  curves~{itc}=L_10;
  If (curves~{itc} < 0)
    curves~{itc}=-L_10;
  EndIf
  itc++;
  If (builtL_11 == -1)
    L_11=newl;
    Line(L_11)={P_12,P_11};
  EndIf
  If (builtL_11 != -1)
    L_11=builtL_11;
  EndIf
  curves~{itc}=L_11;
  If (curves~{itc} < 0)
    curves~{itc}=-L_11;
  EndIf
  itc++;
  If (builtL_12 == -1)
    L_12=newl;
    Line(L_12)={P_11,P_8};
  EndIf
  If (builtL_12 != -1)
    L_12=builtL_12;
  EndIf
  curves~{itc}=L_12;
  If (curves~{itc} < 0)
    curves~{itc}=-L_12;
  EndIf
  itc++;
  If (builtL_13 == -1)
    L_13=newl;
    Line(L_13)={P_8,P_3};
  EndIf
  If (builtL_13 != -1)
    L_13=builtL_13;
  EndIf
  curves~{itc}=L_13;
  If (curves~{itc} < 0)
    curves~{itc}=-L_13;
  EndIf
  itc++;
  If (builtL_14 == -1)
    L_14=newl;
    Line(L_14)={P_10,P_4};
  EndIf
  If (builtL_14 != -1)
    L_14=builtL_14;
  EndIf
  curves~{itc}=L_14;
  If (curves~{itc} < 0)
    curves~{itc}=-L_14;
  EndIf
  itc++;
  If (builtL_15 == -1)
    L_15=newl;
    Line(L_15)={P_1,P_5};
  EndIf
  If (builtL_15 != -1)
    L_15=builtL_15;
  EndIf
  curves~{itc}=L_15;
  If (curves~{itc} < 0)
    curves~{itc}=-L_15;
  EndIf
  itc++;
  If (builtL_16 == -1)
    L_16=newl;
    Line(L_16)={P_2,P_6};
  EndIf
  If (builtL_16 != -1)
    L_16=builtL_16;
  EndIf
  curves~{itc}=L_16;
  If (curves~{itc} < 0)
    curves~{itc}=-L_16;
  EndIf
  itc++;
  If (builtL_17 == -1)
    L_17=newl;
    Line(L_17)={P_7,P_11};
  EndIf
  If (builtL_17 != -1)
    L_17=builtL_17;
  EndIf
  curves~{itc}=L_17;
  If (curves~{itc} < 0)
    curves~{itc}=-L_17;
  EndIf
  itc++;
  If (builtL_18 == -1)
    L_18=newl;
    Line(L_18)={P_9,P_12};
  EndIf
  If (builtL_18 != -1)
    L_18=builtL_18;
  EndIf
  curves~{itc}=L_18;
  If (curves~{itc} < 0)
    curves~{itc}=-L_18;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-L_1,-L_2,-L_3,-L_4};
  LL_2=newll;
  Line Loop(LL_2)={L_5,L_6,L_7,L_8};
  LL_3=newll;
  Line Loop(LL_3)={L_9,L_14,-L_5,-L_13};
  LL_4=newll;
  Line Loop(LL_4)={L_3,L_16,-L_7,-L_15};
  LL_5=newll;
  Line Loop(LL_5)={L_4,L_17,L_12,L_13,-L_8,-L_16};
  LL_6=newll;
  Line Loop(LL_6)={L_2,L_15,-L_6,-L_14,L_10,-L_18};
  LL_7=newll;
  Line Loop(LL_7)={-L_9,-L_10,-L_11,-L_12};
  LL_8=newll;
  Line Loop(LL_8)={L_1,L_18,L_11,-L_17};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8};
Return

Function xlifepp_Parallelepiped5
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;
  If (builtP_9 == -1)
    P_9=newp;
    Point(P_9)={x9,y9,z9,h9};
  EndIf
  If (builtP_9 != -1)
    P_9=builtP_9;
  EndIf
  points~{itp}=P_9;
  itp++;
  If (builtP_10 == -1)
    P_10=newp;
    Point(P_10)={x10,y10,z10,h10};
  EndIf
  If (builtP_10 != -1)
    P_10=builtP_10;
  EndIf
  points~{itp}=P_10;
  itp++;
  If (builtP_11 == -1)
    P_11=newp;
    Point(P_11)={x11,y11,z11,h11};
  EndIf
  If (builtP_11 != -1)
    P_11=builtP_11;
  EndIf
  points~{itp}=P_11;
  itp++;
  If (builtP_12 == -1)
    P_12=newp;
    Point(P_12)={x12,y12,z12,h12};
  EndIf
  If (builtP_12 != -1)
    P_12=builtP_12;
  EndIf
  points~{itp}=P_12;
  itp++;
  If (builtP_13 == -1)
    P_13=newp;
    Point(P_13)={x13,y13,z13,h13};
  EndIf
  If (builtP_13 != -1)
    P_13=builtP_13;
  EndIf
  points~{itp}=P_13;
  itp++;
  If (builtP_14 == -1)
    P_14=newp;
    Point(P_14)={x14,y14,z14,h14};
  EndIf
  If (builtP_14 != -1)
    P_14=builtP_14;
  EndIf
  points~{itp}=P_14;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_11,P_7};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_7,P_1};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_9};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_9,P_11};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_2,P_3};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;
  If (builtL_6 == -1)
    L_6=newl;
    Line(L_6)={P_3,P_4};
  EndIf
  If (builtL_6 != -1)
    L_6=builtL_6;
  EndIf
  curves~{itc}=L_6;
  If (curves~{itc} < 0)
    curves~{itc}=-L_6;
  EndIf
  itc++;
  If (builtL_7 == -1)
    L_7=newl;
    Line(L_7)={P_4,P_5};
  EndIf
  If (builtL_7 != -1)
    L_7=builtL_7;
  EndIf
  curves~{itc}=L_7;
  If (curves~{itc} < 0)
    curves~{itc}=-L_7;
  EndIf
  itc++;
  If (builtL_8 == -1)
    L_8=newl;
    Line(L_8)={P_5,P_2};
  EndIf
  If (builtL_8 != -1)
    L_8=builtL_8;
  EndIf
  curves~{itc}=L_8;
  If (curves~{itc} < 0)
    curves~{itc}=-L_8;
  EndIf
  itc++;
  If (builtL_9 == -1)
    L_9=newl;
    Line(L_9)={P_6,P_8};
  EndIf
  If (builtL_9 != -1)
    L_9=builtL_9;
  EndIf
  curves~{itc}=L_9;
  If (curves~{itc} < 0)
    curves~{itc}=-L_9;
  EndIf
  itc++;
  If (builtL_10 == -1)
    L_10=newl;
    Line(L_10)={P_8,P_12};
  EndIf
  If (builtL_10 != -1)
    L_10=builtL_10;
  EndIf
  curves~{itc}=L_10;
  If (curves~{itc} < 0)
    curves~{itc}=-L_10;
  EndIf
  itc++;
  If (builtL_11 == -1)
    L_11=newl;
    Line(L_11)={P_12,P_14};
  EndIf
  If (builtL_11 != -1)
    L_11=builtL_11;
  EndIf
  curves~{itc}=L_11;
  If (curves~{itc} < 0)
    curves~{itc}=-L_11;
  EndIf
  itc++;
  If (builtL_12 == -1)
    L_12=newl;
    Line(L_12)={P_14,P_13};
  EndIf
  If (builtL_12 != -1)
    L_12=builtL_12;
  EndIf
  curves~{itc}=L_12;
  If (curves~{itc} < 0)
    curves~{itc}=-L_12;
  EndIf
  itc++;
  If (builtL_13 == -1)
    L_13=newl;
    Line(L_13)={P_13,P_10};
  EndIf
  If (builtL_13 != -1)
    L_13=builtL_13;
  EndIf
  curves~{itc}=L_13;
  If (curves~{itc} < 0)
    curves~{itc}=-L_13;
  EndIf
  itc++;
  If (builtL_14 == -1)
    L_14=newl;
    Line(L_14)={P_10,P_6};
  EndIf
  If (builtL_14 != -1)
    L_14=builtL_14;
  EndIf
  curves~{itc}=L_14;
  If (curves~{itc} < 0)
    curves~{itc}=-L_14;
  EndIf
  itc++;
  If (builtL_15 == -1)
    L_15=newl;
    Line(L_15)={P_6,P_2};
  EndIf
  If (builtL_15 != -1)
    L_15=builtL_15;
  EndIf
  curves~{itc}=L_15;
  If (curves~{itc} < 0)
    curves~{itc}=-L_15;
  EndIf
  itc++;
  If (builtL_16 == -1)
    L_16=newl;
    Line(L_16)={P_8,P_3};
  EndIf
  If (builtL_16 != -1)
    L_16=builtL_16;
  EndIf
  curves~{itc}=L_16;
  If (curves~{itc} < 0)
    curves~{itc}=-L_16;
  EndIf
  itc++;
  If (builtL_17 == -1)
    L_17=newl;
    Line(L_17)={P_1,P_4};
  EndIf
  If (builtL_17 != -1)
    L_17=builtL_17;
  EndIf
  curves~{itc}=L_17;
  If (curves~{itc} < 0)
    curves~{itc}=-L_17;
  EndIf
  itc++;
  If (builtL_18 == -1)
    L_18=newl;
    Line(L_18)={P_10,P_5};
  EndIf
  If (builtL_18 != -1)
    L_18=builtL_18;
  EndIf
  curves~{itc}=L_18;
  If (curves~{itc} < 0)
    curves~{itc}=-L_18;
  EndIf
  itc++;
  If (builtL_19 == -1)
    L_19=newl;
    Line(L_19)={P_11,P_14};
  EndIf
  If (builtL_19 != -1)
    L_19=builtL_19;
  EndIf
  curves~{itc}=L_19;
  If (curves~{itc} < 0)
    curves~{itc}=-L_19;
  EndIf
  itc++;
  If (builtL_20 == -1)
    L_20=newl;
    Line(L_20)={P_7,P_12};
  EndIf
  If (builtL_20 != -1)
    L_20=builtL_20;
  EndIf
  curves~{itc}=L_20;
  If (curves~{itc} < 0)
    curves~{itc}=-L_20;
  EndIf
  itc++;
  If (builtL_21 == -1)
    L_21=newl;
    Line(L_21)={P_9,P_13};
  EndIf
  If (builtL_21 != -1)
    L_21=builtL_21;
  EndIf
  curves~{itc}=L_21;
  If (curves~{itc} < 0)
    curves~{itc}=-L_21;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-L_1,-L_2,-L_3,-L_4};
  LL_2=newll;
  Line Loop(LL_2)={L_5,L_6,L_7,L_8};
  LL_3=newll;
  Line Loop(LL_3)={L_9,L_16,-L_5,-L_15};
  LL_4=newll;
  Line Loop(LL_4)={L_3,L_21,L_13,L_18,-L_7,-L_17};
  LL_5=newll;
  Line Loop(LL_5)={L_14,L_15,-L_8,-L_18};
  LL_6=newll;
  Line Loop(LL_6)={L_2,L_17,-L_6,-L_16,L_10,-L_20};
  LL_7=newll;
  Line Loop(LL_7)={-L_9,-L_14,-L_13,-L_12,-L_11,-L_10};
  LL_8=newll;
  Line Loop(LL_8)={L_1,L_20,L_11,-L_19};
  LL_9=newll;
  Line Loop(LL_9)={L_4,L_19,L_12,-L_21};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;
  If (builtS_9 == -1)
    S_9=news;
    Plane Surface(S_9)={LL_9};
  EndIf
  If (builtS_9 != -1)
    S_9=builtS_9;
  EndIf
  surfs~{its}=S_9;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8,S_9};
Return

Function xlifepp_Parallelepiped3
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  points~{itp}=P_1;
  itp++;
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;
  If (builtP_9 == -1)
    P_9=newp;
    Point(P_9)={x9,y9,z9,h9};
  EndIf
  If (builtP_9 != -1)
    P_9=builtP_9;
  EndIf
  points~{itp}=P_9;
  itp++;
  If (builtP_10 == -1)
    P_10=newp;
    Point(P_10)={x10,y10,z10,h10};
  EndIf
  If (builtP_10 != -1)
    P_10=builtP_10;
  EndIf
  points~{itp}=P_10;
  itp++;
  If (builtP_11 == -1)
    P_11=newp;
    Point(P_11)={x11,y11,z11,h11};
  EndIf
  If (builtP_11 != -1)
    P_11=builtP_11;
  EndIf
  points~{itp}=P_11;
  itp++;
  If (builtP_12 == -1)
    P_12=newp;
    Point(P_12)={x12,y12,z12,h12};
  EndIf
  If (builtP_12 != -1)
    P_12=builtP_12;
  EndIf
  points~{itp}=P_12;
  itp++;

  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_4,P_9};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_9,P_12};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_12,P_10};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_10,P_5};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_5,P_6};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;
  If (builtL_6 == -1)
    L_6=newl;
    Line(L_6)={P_6,P_4};
  EndIf
  If (builtL_6 != -1)
    L_6=builtL_6;
  EndIf
  curves~{itc}=L_6;
  If (curves~{itc} < 0)
    curves~{itc}=-L_6;
  EndIf
  itc++;
  If (builtL_7 == -1)
    L_7=newl;
    Line(L_7)={P_1,P_7};
  EndIf
  If (builtL_7 != -1)
    L_7=builtL_7;
  EndIf
  curves~{itc}=L_7;
  If (curves~{itc} < 0)
    curves~{itc}=-L_7;
  EndIf
  itc++;
  If (builtL_8 == -1)
    L_8=newl;
    Line(L_8)={P_7,P_11};
  EndIf
  If (builtL_8 != -1)
    L_8=builtL_8;
  EndIf
  curves~{itc}=L_8;
  If (curves~{itc} < 0)
    curves~{itc}=-L_8;
  EndIf
  itc++;
  If (builtL_9 == -1)
    L_9=newl;
    Line(L_9)={P_11,P_8};
  EndIf
  If (builtL_9 != -1)
    L_9=builtL_9;
  EndIf
  curves~{itc}=L_9;
  If (curves~{itc} < 0)
    curves~{itc}=-L_9;
  EndIf
  itc++;
  If (builtL_10 == -1)
    L_10=newl;
    Line(L_10)={P_8,P_2};
  EndIf
  If (builtL_10 != -1)
    L_10=builtL_10;
  EndIf
  curves~{itc}=L_10;
  If (curves~{itc} < 0)
    curves~{itc}=-L_10;
  EndIf
  itc++;
  If (builtL_11 == -1)
    L_11=newl;
    Line(L_11)={P_2,P_3};
  EndIf
  If (builtL_11 != -1)
    L_11=builtL_11;
  EndIf
  curves~{itc}=L_11;
  If (curves~{itc} < 0)
    curves~{itc}=-L_11;
  EndIf
  itc++;
  If (builtL_12 == -1)
    L_12=newl;
    Line(L_12)={P_3,P_1};
  EndIf
  If (builtL_12 != -1)
    L_12=builtL_12;
  EndIf
  curves~{itc}=L_12;
  If (curves~{itc} < 0)
    curves~{itc}=-L_12;
  EndIf
  itc++;
  If (builtL_13 == -1)
    L_13=newl;
    Line(L_13)={P_4,P_1};
  EndIf
  If (builtL_13 != -1)
    L_13=builtL_13;
  EndIf
  curves~{itc}=L_13;
  If (curves~{itc} < 0)
    curves~{itc}=-L_13;
  EndIf
  itc++;
  If (builtL_14 == -1)
    L_14=newl;
    Line(L_14)={P_9,P_7};
  EndIf
  If (builtL_14 != -1)
    L_14=builtL_14;
  EndIf
  curves~{itc}=L_14;
  If (curves~{itc} < 0)
    curves~{itc}=-L_14;
  EndIf
  itc++;
  If (builtL_15 == -1)
    L_15=newl;
    Line(L_15)={P_12,P_11};
  EndIf
  If (builtL_15 != -1)
    L_15=builtL_15;
  EndIf
  curves~{itc}=L_15;
  If (curves~{itc} < 0)
    curves~{itc}=-L_15;
  EndIf
  itc++;
  If (builtL_16 == -1)
    L_16=newl;
    Line(L_16)={P_10,P_8};
  EndIf
  If (builtL_16 != -1)
    L_16=builtL_16;
  EndIf
  curves~{itc}=L_16;
  If (curves~{itc} < 0)
    curves~{itc}=-L_16;
  EndIf
  itc++;
  If (builtL_17 == -1)
    L_17=newl;
    Line(L_17)={P_5,P_2};
  EndIf
  If (builtL_17 != -1)
    L_17=builtL_17;
  EndIf
  curves~{itc}=L_17;
  If (curves~{itc} < 0)
    curves~{itc}=-L_17;
  EndIf
  itc++;
  If (builtL_18 == -1)
    L_18=newl;
    Line(L_18)={P_6,P_3};
  EndIf
  If (builtL_18 != -1)
    L_18=builtL_18;
  EndIf
  curves~{itc}=L_18;
  If (curves~{itc} < 0)
    curves~{itc}=-L_18;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-L_1,-L_2,-L_3,-L_4,-L_5,-L_6};
  LL_2=newll;
  Line Loop(LL_2)={L_7,L_8,L_9,L_10,L_11,L_12};
  LL_3=newll;
  Line Loop(LL_3)={L_1,L_14,-L_7,-L_13};
  LL_4=newll;
  Line Loop(LL_4)={L_5,L_18,-L_11,-L_17};
  LL_5=newll;
  Line Loop(LL_5)={L_6,L_13,-L_12,-L_18};
  LL_6=newll;
  Line Loop(LL_6)={L_4,L_17,-L_10,-L_16};
  LL_7=newll;
  Line Loop(LL_7)={L_3,L_16,-L_9,-L_15};
  LL_8=newll;
  Line Loop(LL_8)={L_2,L_15,-L_8,-L_14};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8};
Return

Function xlifepp_Ellipsoid
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_2,P_1,P_2,P_7};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_7,P_1,P_4,P_4};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_4,P_1,P_4,P_6};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_6,P_1,P_2,P_2};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtE_9 == -1)
    E_9=newl;
    Ellipse(E_9)={P_3,P_1,P_3,P_7};
  EndIf
  If (builtE_9 != -1)
    E_9=builtE_9;
  EndIf
  curves~{itc}=E_9;
  If (curves~{itc} < 0)
    curves~{itc}=-E_9;
  EndIf
  itc++;
  If (builtE_10 == -1)
    E_10=newl;
    Ellipse(E_10)={P_7,P_1,P_5,P_5};
  EndIf
  If (builtE_10 != -1)
    E_10=builtE_10;
  EndIf
  curves~{itc}=E_10;
  If (curves~{itc} < 0)
    curves~{itc}=-E_10;
  EndIf
  itc++;
  If (builtE_11 == -1)
    E_11=newl;
    Ellipse(E_11)={P_5,P_1,P_5,P_6};
  EndIf
  If (builtE_11 != -1)
    E_11=builtE_11;
  EndIf
  curves~{itc}=E_11;
  If (curves~{itc} < 0)
    curves~{itc}=-E_11;
  EndIf
  itc++;
  If (builtE_12 == -1)
    E_12=newl;
    Ellipse(E_12)={P_6,P_1,P_3,P_3};
  EndIf
  If (builtE_12 != -1)
    E_12=builtE_12;
  EndIf
  curves~{itc}=E_12;
  If (curves~{itc} < 0)
    curves~{itc}=-E_12;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_9,-E_5};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_9,-E_6};
  LL_3=newll;
  Line Loop(LL_3)={E_3,-E_10,E_6};
  LL_4=newll;
  Line Loop(LL_4)={E_4,E_10,E_5};
  LL_5=newll;
  Line Loop(LL_5)={-E_1,E_12,-E_8};
  LL_6=newll;
  Line Loop(LL_6)={-E_2,-E_12,-E_7};
  LL_7=newll;
  Line Loop(LL_7)={-E_3,-E_11,E_7};
  LL_8=newll;
  Line Loop(LL_8)={-E_4,E_11,E_8};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Ruled Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Ruled Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Ruled Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8};
Return

Function xlifepp_Ellipsoid7
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_2,P_1,P_2,P_7};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_7,P_1,P_4,P_4};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_4,P_1,P_4,P_6};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_6,P_1,P_2,P_2};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtE_9 == -1)
    E_9=newl;
    Ellipse(E_9)={P_3,P_1,P_3,P_7};
  EndIf
  If (builtE_9 != -1)
    E_9=builtE_9;
  EndIf
  curves~{itc}=E_9;
  If (curves~{itc} < 0)
    curves~{itc}=-E_9;
  EndIf
  itc++;
  If (builtE_10 == -1)
    E_10=newl;
    Ellipse(E_10)={P_7,P_1,P_5,P_5};
  EndIf
  If (builtE_10 != -1)
    E_10=builtE_10;
  EndIf
  curves~{itc}=E_10;
  If (curves~{itc} < 0)
    curves~{itc}=-E_10;
  EndIf
  itc++;
  If (builtE_11 == -1)
    E_11=newl;
    Ellipse(E_11)={P_5,P_1,P_5,P_6};
  EndIf
  If (builtE_11 != -1)
    E_11=builtE_11;
  EndIf
  curves~{itc}=E_11;
  If (curves~{itc} < 0)
    curves~{itc}=-E_11;
  EndIf
  itc++;
  If (builtE_12 == -1)
    E_12=newl;
    Ellipse(E_12)={P_6,P_1,P_3,P_3};
  EndIf
  If (builtE_12 != -1)
    E_12=builtE_12;
  EndIf
  curves~{itc}=E_12;
  If (curves~{itc} < 0)
    curves~{itc}=-E_12;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_5};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_6};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_9,-E_5};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_9,-E_6};
  LL_3=newll;
  Line Loop(LL_3)={E_3,-E_10,E_6};
  LL_4=newll;
  Line Loop(LL_4)={E_4,E_10,E_5};
  LL_5=newll;
  Line Loop(LL_5)={-E_1,E_12,-E_8};
  LL_6=newll;
  Line Loop(LL_6)={-E_2,-E_12,-E_7};
  LL_7=newll;
  Line Loop(LL_7)={-E_3,-E_11,E_7};
  LL_8=newll;
  Line Loop(LL_8)={L_2,E_11,-L_3};
  LL_9=newll;
  Line Loop(LL_9)={L_3,E_8,-L_1};
  LL_10=newll;
  Line Loop(LL_10)={L_1,-E_4,-L_2};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Ruled Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Ruled Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;
  If (builtS_9 == -1)
    S_9=news;
    Plane Surface(S_9)={LL_9};
  EndIf
  If (builtS_9 != -1)
    S_9=builtS_9;
  EndIf
  surfs~{its}=S_9;
  its++;
  If (builtS_10 == -1)
    S_10=news;
    Plane Surface(S_10)={LL_10};
  EndIf
  If (builtS_10 != -1)
    S_10=builtS_10;
  EndIf
  surfs~{its}=S_10;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8,S_9,S_10};
Return

Function xlifepp_Ellipsoid6
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_2,P_1,P_2,P_7};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_7,P_1,P_4,P_4};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_4,P_1,P_4,P_6};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_6,P_1,P_2,P_2};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtE_9 == -1)
    E_9=newl;
    Ellipse(E_9)={P_3,P_1,P_3,P_7};
  EndIf
  If (builtE_9 != -1)
    E_9=builtE_9;
  EndIf
  curves~{itc}=E_9;
  If (curves~{itc} < 0)
    curves~{itc}=-E_9;
  EndIf
  itc++;
  If (builtE_10 == -1)
    E_10=newl;
    Ellipse(E_10)={P_7,P_1,P_5,P_5};
  EndIf
  If (builtE_10 != -1)
    E_10=builtE_10;
  EndIf
  curves~{itc}=E_10;
  If (curves~{itc} < 0)
    curves~{itc}=-E_10;
  EndIf
  itc++;
  If (builtE_11 == -1)
    E_11=newl;
    Ellipse(E_11)={P_6,P_1,P_3,P_3};
  EndIf
  If (builtE_11 != -1)
    E_11=builtE_11;
  EndIf
  curves~{itc}=E_11;
  If (curves~{itc} < 0)
    curves~{itc}=-E_11;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_4};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_5};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_1,P_6};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_9,-E_5};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_9,-E_6};
  LL_3=newll;
  Line Loop(LL_3)={E_3,-E_10,E_6};
  LL_4=newll;
  Line Loop(LL_4)={E_4,E_10,E_5};
  LL_5=newll;
  Line Loop(LL_5)={-E_1,E_11,-E_8};
  LL_6=newll;
  Line Loop(LL_6)={-E_2,-E_11,-E_7};
  LL_7=newll;
  Line Loop(LL_7)={L_2,E_7,-L_4};
  LL_8=newll;
  Line Loop(LL_8)={L_4,E_8,-L_1};
  LL_9=newll;
  Line Loop(LL_9)={L_3,-E_3,-L_2};
  LL_10=newll;
  Line Loop(LL_10)={L_1,-E_4,-L_3};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Ruled Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;
  If (builtS_9 == -1)
    S_9=news;
    Plane Surface(S_9)={LL_9};
  EndIf
  If (builtS_9 != -1)
    S_9=builtS_9;
  EndIf
  surfs~{its}=S_9;
  its++;
  If (builtS_10 == -1)
    S_10=news;
    Plane Surface(S_10)={LL_10};
  EndIf
  If (builtS_10 != -1)
    S_10=builtS_10;
  EndIf
  surfs~{its}=S_10;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8,S_9,S_10};
Return

Function xlifepp_Ellipsoid5
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_2,P_1,P_2,P_7};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_7,P_1,P_4,P_4};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_6,P_1,P_2,P_2};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_3,P_1,P_3,P_7};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtE_9 == -1)
    E_9=newl;
    Ellipse(E_9)={P_7,P_1,P_5,P_5};
  EndIf
  If (builtE_9 != -1)
    E_9=builtE_9;
  EndIf
  curves~{itc}=E_9;
  If (curves~{itc} < 0)
    curves~{itc}=-E_9;
  EndIf
  itc++;
  If (builtE_10 == -1)
    E_10=newl;
    Ellipse(E_10)={P_6,P_1,P_3,P_3};
  EndIf
  If (builtE_10 != -1)
    E_10=builtE_10;
  EndIf
  curves~{itc}=E_10;
  If (curves~{itc} < 0)
    curves~{itc}=-E_10;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_4};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_3};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_1,P_5};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_1,P_6};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_8,-E_5};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_8,-E_6};
  LL_3=newll;
  Line Loop(LL_3)={E_3,-E_9,E_6};
  LL_4=newll;
  Line Loop(LL_4)={E_4,E_9,E_5};
  LL_5=newll;
  Line Loop(LL_5)={-E_1,E_10,-E_7};
  LL_6=newll;
  Line Loop(LL_6)={L_3,-E_10,-L_5};
  LL_7=newll;
  Line Loop(LL_7)={L_5,E_7,-L_1};
  LL_8=newll;
  Line Loop(LL_8)={L_2,-E_2,-L_3};
  LL_9=newll;
  Line Loop(LL_9)={L_4,-E_3,-L_2};
  LL_10=newll;
  Line Loop(LL_10)={L_1,-E_4,-L_4};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;
  If (builtS_9 == -1)
    S_9=news;
    Plane Surface(S_9)={LL_9};
  EndIf
  If (builtS_9 != -1)
    S_9=builtS_9;
  EndIf
  surfs~{its}=S_9;
  its++;
  If (builtS_10 == -1)
    S_10=news;
    Plane Surface(S_10)={LL_10};
  EndIf
  If (builtS_10 != -1)
    S_10=builtS_10;
  EndIf
  surfs~{its}=S_10;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8,S_9,S_10};
Return

Function xlifepp_Ellipsoid4
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_2,P_1,P_2,P_6};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_6,P_1,P_4,P_4};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_3,P_1,P_3,P_6};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_6,P_1,P_5,P_5};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_9 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_4};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_3};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_1,P_5};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_7,-E_5};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_7,-E_6};
  LL_3=newll;
  Line Loop(LL_3)={E_3,-E_8,E_6};
  LL_4=newll;
  Line Loop(LL_4)={E_4,E_8,E_5};
  LL_5=newll;
  Line Loop(LL_5)={L_3,-E_1,-L_1};
  LL_6=newll;
  Line Loop(LL_6)={L_2,-E_2,-L_3};
  LL_7=newll;
  Line Loop(LL_7)={L_4,-E_3,-L_2};
  LL_8=newll;
  Line Loop(LL_8)={L_1,-E_4,-L_4};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Plane Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Plane Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8};
Return

Function xlifepp_Ellipsoid3
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  points~{itp}=P_6;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_2,P_1,P_2,P_6};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_6,P_1,P_4,P_4};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_3,P_1,P_3,P_6};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_6,P_1,P_5,P_5};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_4};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_3};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_1,P_5};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;
  If (builtL_5 == -1)
    L_5=newl;
    Line(L_5)={P_1,P_6};
  EndIf
  If (builtL_5 != -1)
    L_5=builtL_5;
  EndIf
  curves~{itc}=L_5;
  If (curves~{itc} < 0)
    curves~{itc}=-L_5;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_6,-E_4};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_6,-E_5};
  LL_3=newll;
  Line Loop(LL_3)={E_3,-E_7,E_5};
  LL_4=newll;
  Line Loop(LL_4)={E_7,-L_4,L_5};
  LL_5=newll;
  Line Loop(LL_5)={L_1,E_4,-L_5};
  LL_6=newll;
  Line Loop(LL_6)={L_3,-E_1,-L_1};
  LL_7=newll;
  Line Loop(LL_7)={L_2,-E_2,-L_3};
  LL_8=newll;
  Line Loop(LL_8)={L_4,-E_3,-L_2};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Ruled Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;
  If (builtS_7 == -1)
    S_7=news;
    Ruled Surface(S_7)={LL_7};
  EndIf
  If (builtS_7 != -1)
    S_7=builtS_7;
  EndIf
  surfs~{its}=S_7;
  its++;
  If (builtS_8 == -1)
    S_8=news;
    Ruled Surface(S_8)={LL_8};
  EndIf
  If (builtS_8 != -1)
    S_8=builtS_8;
  EndIf
  surfs~{its}=S_8;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8};
Return

Function xlifepp_Ellipsoid2
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_2,P_1,P_2,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_4,P_4};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_3,P_1,P_5,P_5};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_4};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_3};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_1,P_5};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_5,-E_3};
  LL_2=newll;
  Line Loop(LL_2)={E_2,-E_5,-E_4};
  LL_3=newll;
  Line Loop(LL_3)={E_4,-L_2,L_4};
  LL_4=newll;
  Line Loop(LL_4)={E_3,-L_4,L_1};
  LL_5=newll;
  Line Loop(LL_5)={L_3,-E_1,-L_1};
  LL_6=newll;
  Line Loop(LL_6)={L_2,-E_2,-L_3};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Plane Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Plane Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Plane Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Plane Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6};
Return

Function xlifepp_Ellipsoid1
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_2,P_1,P_2,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 == -1)
    E_3=newl;
    Ellipse(E_3)={P_3,P_1,P_3,P_4};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_1,P_2};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_1,P_3};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_1,P_4};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={E_1,E_3,-E_2};
  LL_2=newll;
  Line Loop(LL_2)={L_3,-E_3,-L_2};
  LL_3=newll;
  Line Loop(LL_3)={L_1,E_2,-L_3};
  LL_4=newll;
  Line Loop(LL_4)={L_2,-E_1,-L_1};

  If (builtS_1 == -1)
    S_1=news;
    Ruled Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Ruled Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If (builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4};
Return

Function xlifepp_RevTrunk
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;
  If (builtP_9 == -1)
    P_9=newp;
    Point(P_9)={x9,y9,z9,h9};
  EndIf
  If (builtP_9 != -1)
    P_9=builtP_9;
  EndIf
  points~{itp}=P_9;
  itp++;
  If (builtP_10 == -1)
    P_10=newp;
    Point(P_10)={x10,y10,z10,h10};
  EndIf
  If (builtP_10 != -1)
    P_10=builtP_10;
  EndIf
  points~{itp}=P_10;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 ==-1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_7,P_6,P_7,P_8};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_8,P_6,P_9,P_9};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_9,P_6,P_9,P_10};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_10,P_6,P_7,P_7};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_2,P_7};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_3,P_8};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_4,P_9};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_5,P_10};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-E_1,-E_2,-E_3,-E_4};
  LL_2=newll;
  Line Loop(LL_2)={E_5,E_6,E_7,E_8};
  LL_3=newll;
  Line Loop(LL_3)={E_1,L_2,-E_5,-L_1};
  LL_4=newll;
  Line Loop(LL_4)={E_2,L_3,-E_6,-L_2};
  LL_5=newll;
  Line Loop(LL_5)={E_3,L_4,-E_7,-L_3};
  LL_6=newll;
  Line Loop(LL_6)={E_4,L_1,-E_8,-L_4};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If(builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Ruled Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6};
Return

Function xlifepp_RevCylinder
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  points~{itp}=P_2;
  itp++;
  If (builtP_3 == -1)
    P_3=newp;
    Point(P_3)={x3,y3,z3,h3};
  EndIf
  If (builtP_3 != -1)
    P_3=builtP_3;
  EndIf
  points~{itp}=P_3;
  itp++;
  If (builtP_4 == -1)
    P_4=newp;
    Point(P_4)={x4,y4,z4,h4};
  EndIf
  If (builtP_4 != -1)
    P_4=builtP_4;
  EndIf
  points~{itp}=P_4;
  itp++;
  If (builtP_5 == -1)
    P_5=newp;
    Point(P_5)={x5,y5,z5,h5};
  EndIf
  If (builtP_5 != -1)
    P_5=builtP_5;
  EndIf
  points~{itp}=P_5;
  itp++;
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf
  If (builtP_6 != -1)
    P_6=builtP_6;
  EndIf
  If (builtP_7 == -1)
    P_7=newp;
    Point(P_7)={x7,y7,z7,h7};
  EndIf
  If (builtP_7 != -1)
    P_7=builtP_7;
  EndIf
  points~{itp}=P_7;
  itp++;
  If (builtP_8 == -1)
    P_8=newp;
    Point(P_8)={x8,y8,z8,h8};
  EndIf
  If (builtP_8 != -1)
    P_8=builtP_8;
  EndIf
  points~{itp}=P_8;
  itp++;
  If (builtP_9 == -1)
    P_9=newp;
    Point(P_9)={x9,y9,z9,h9};
  EndIf
  If (builtP_9 != -1)
    P_9=builtP_9;
  EndIf
  points~{itp}=P_9;
  itp++;
  If (builtP_10 == -1)
    P_10=newp;
    Point(P_10)={x10,y10,z10,h10};
  EndIf
  If (builtP_10 != -1)
    P_10=builtP_10;
  EndIf
  points~{itp}=P_10;
  itp++;

  If (builtE_1 == -1)
    E_1=newl;
    Ellipse(E_1)={P_2,P_1,P_2,P_3};
  EndIf
  If (builtE_1 != -1)
    E_1=builtE_1;
  EndIf
  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  If (builtE_2 == -1)
    E_2=newl;
    Ellipse(E_2)={P_3,P_1,P_4,P_4};
  EndIf
  If (builtE_2 != -1)
    E_2=builtE_2;
  EndIf
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  If (builtE_3 ==-1)
    E_3=newl;
    Ellipse(E_3)={P_4,P_1,P_4,P_5};
  EndIf
  If (builtE_3 != -1)
    E_3=builtE_3;
  EndIf
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  If (builtE_4 == -1)
    E_4=newl;
    Ellipse(E_4)={P_5,P_1,P_2,P_2};
  EndIf
  If (builtE_4 != -1)
    E_4=builtE_4;
  EndIf
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  If (builtE_5 == -1)
    E_5=newl;
    Ellipse(E_5)={P_7,P_6,P_7,P_8};
  EndIf
  If (builtE_5 != -1)
    E_5=builtE_5;
  EndIf
  curves~{itc}=E_5;
  If (curves~{itc} < 0)
    curves~{itc}=-E_5;
  EndIf
  itc++;
  If (builtE_6 == -1)
    E_6=newl;
    Ellipse(E_6)={P_8,P_6,P_9,P_9};
  EndIf
  If (builtE_6 != -1)
    E_6=builtE_6;
  EndIf
  curves~{itc}=E_6;
  If (curves~{itc} < 0)
    curves~{itc}=-E_6;
  EndIf
  itc++;
  If (builtE_7 == -1)
    E_7=newl;
    Ellipse(E_7)={P_9,P_6,P_9,P_10};
  EndIf
  If (builtE_7 != -1)
    E_7=builtE_7;
  EndIf
  curves~{itc}=E_7;
  If (curves~{itc} < 0)
    curves~{itc}=-E_7;
  EndIf
  itc++;
  If (builtE_8 == -1)
    E_8=newl;
    Ellipse(E_8)={P_10,P_6,P_7,P_7};
  EndIf
  If (builtE_8 != -1)
    E_8=builtE_8;
  EndIf
  curves~{itc}=E_8;
  If (curves~{itc} < 0)
    curves~{itc}=-E_8;
  EndIf
  itc++;
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_2,P_7};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  If (builtL_2 == -1)
    L_2=newl;
    Line(L_2)={P_3,P_8};
  EndIf
  If (builtL_2 != -1)
    L_2=builtL_2;
  EndIf
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  If (builtL_3 == -1)
    L_3=newl;
    Line(L_3)={P_4,P_9};
  EndIf
  If (builtL_3 != -1)
    L_3=builtL_3;
  EndIf
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  If (builtL_4 == -1)
    L_4=newl;
    Line(L_4)={P_5,P_10};
  EndIf
  If (builtL_4 != -1)
    L_4=builtL_4;
  EndIf
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  LL_1=newll;
  Line Loop(LL_1)={-E_1,-E_2,-E_3,-E_4};
  LL_2=newll;
  Line Loop(LL_2)={E_5,E_6,E_7,E_8};
  LL_3=newll;
  Line Loop(LL_3)={E_1,L_2,-E_5,-L_1};
  LL_4=newll;
  Line Loop(LL_4)={E_2,L_3,-E_6,-L_2};
  LL_5=newll;
  Line Loop(LL_5)={E_3,L_4,-E_7,-L_3};
  LL_6=newll;
  Line Loop(LL_6)={E_4,L_1,-E_8,-L_4};

  If (builtS_1 == -1)
    S_1=news;
    Plane Surface(S_1)={LL_1};
  EndIf
  If (builtS_1 != -1)
    S_1=builtS_1;
  EndIf
  surfs~{its}=S_1;
  its++;
  If (builtS_2 == -1)
    S_2=news;
    Plane Surface(S_2)={LL_2};
  EndIf
  If (builtS_2 != -1)
    S_2=builtS_2;
  EndIf
  surfs~{its}=S_2;
  its++;
  If (builtS_3 == -1)
    S_3=news;
    Ruled Surface(S_3)={LL_3};
  EndIf
  If (builtS_3 != -1)
    S_3=builtS_3;
  EndIf
  surfs~{its}=S_3;
  its++;
  If(builtS_4 == -1)
    S_4=news;
    Ruled Surface(S_4)={LL_4};
  EndIf
  If (builtS_4 != -1)
    S_4=builtS_4;
  EndIf
  surfs~{its}=S_4;
  its++;
  If (builtS_5 == -1)
    S_5=news;
    Ruled Surface(S_5)={LL_5};
  EndIf
  If (builtS_5 != -1)
    S_5=builtS_5;
  EndIf
  surfs~{its}=S_5;
  its++;
  If (builtS_6 == -1)
    S_6=news;
    Ruled Surface(S_6)={LL_6};
  EndIf
  If (builtS_6 != -1)
    S_6=builtS_6;
  EndIf
  surfs~{its}=S_6;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5,S_6};
Return

Function xlifepp_RevCone
  If (builtP_1 == -1)
    P_1=newp;
    Point(P_1)={x1,y1,z1,h1};
  EndIf
  If (builtP_1 != -1)
    P_1=builtP_1;
  EndIf
  If (builtP_2 == -1)
    P_2=newp;
    Point(P_2)={x2,y2,z2,h2};
  EndIf
  If (builtP_2 != -1)
    P_2=builtP_2;
  EndIf
  If (builtP_6 == -1)
    P_6=newp;
    Point(P_6)={x6,y6,z6,h6};
  EndIf

  TL_1=newl;
  Line(TL_1)={P_1 ,P_2};
  If (builtL_1 == -1)
    L_1=newl;
    Line(L_1)={P_2 ,P_6};
  EndIf
  If (builtL_1 != -1)
    L_1=builtL_1;
  EndIf
  A1=newl;
  Line(A1)={P_6 ,P_1};

  TLL_1=newll;
  Line Loop(TLL_1)={TL_1,L_1,A1};
  TS_1=news;
  Plane Surface(TS_1)={TLL_1};

  out[]=Extrude { {ux,uy,uz}, {x1,y1,z1}, Pi/2} {Surface{TS_1}; Line{L_1, TL_1}; Point{P_2};};
  TS_2=out[0]; // image of TS_1
  V11=out[1]; // first quarter of volume
  S_11=out[2]; // first quarter of base surface (S_1)
  S_2=out[3]; // first quarter of lateral surface
  L_2=out[4]; // image of L_1
  TL_2=out[5]; // new segment inside the base
  P_3=out[6]; // image of P_2

  // We don't have the number of the circle arc created !!!
  // We use the Boundary command and compare with TL_1 and TL_2 to get the third boundary number
  out2[]=Boundary{Surface{S_11};};
  For i In {0:#out2[]-1}
    If (Fabs(out2[i]) != TL_1)
      If (Fabs(out2[i]) != TL_2)
        E_1=out2[i];
      EndIf
    EndIf
  EndFor

  out[]=Extrude { {ux,uy,uz}, {x1,y1,z1}, Pi/2} {Surface{TS_2}; Line{L_2, TL_2}; Point{P_3};};
  TS_3=out[0]; // image of TS_2
  V12=out[1]; // second quarter of volume
  S_12=out[2]; // second quarter of base surface (S_1)
  S_3=out[3]; // second quarter of lateral surface
  L_3=out[4]; // image of L_2
  TL_3=out[5]; // new segment inside the base
  P_4=out[6]; // image of P_3

  // We don't have the number of the circle arc created !!!
  // We use the Boundary command and compare with TL_2 and TL_3 to get the third boundary number
  out2[]=Boundary{Surface{S_12};};
  For i In {0:#out2[]-1}
    If (Fabs(out2[i]) != TL_2)
      If (Fabs(out2[i]) != TL_3)
        E_2=out2[i];
      EndIf
    EndIf
  EndFor

  out[]=Extrude { {ux,uy,uz}, {x1,y1,z1}, Pi/2} {Surface{TS_3}; Line{L_3, TL_3}; Point{P_4};};
  TS_4=out[0]; // image of TS_3
  V13=out[1]; // third quarter of volume
  S_13=out[2]; // third quarter of base surface (S_1)
  S_4=out[3]; // third quarter of lateral surface
  L_4=out[4]; // image of L_3
  TL_4=out[5]; // new segment inside the base
  P_5=out[6]; // image of P_4

  // We don't have the number of the circle arc created !!!
  // We use the Boundary command and compare with TL_3 and TL_4 to get the third boundary number
  out2[]=Boundary{Surface{S_13};};
  For i In {0:#out2[]-1}
    If (Fabs(out2[i]) != TL_3)
      If (Fabs(out2[i]) != TL_4)
        E_3=out2[i];
      EndIf
    EndIf
  EndFor

  out[]=Extrude { {ux,uy,uz}, {x1,y1,z1}, Pi/2} {Surface{TS_4}; Line{L_4,TL_4}; Point{P_5};};
  If (TS_1 != out[0])
    Printf("ERROR Surface");
  EndIf
  V14=out[1]; // last quarter of volume
  S_14=out[2]; // last quarter of base surface (S_1)
  S_5=out[3]; // last quarter of lateral surface

  // We don't have the number of the circle arc created !!!
  // We use the Boundary command and compare with TL_4 and TL_1 to get the third boundary number
  out2[]=Boundary{Surface{S_14};};
  For i In {0:#out2[]-1}
    If (Fabs(out2[i]) != TL_4)
      If (Fabs(out2[i]) != TL_1)
        E_4=out2[i];
      EndIf
    EndIf
  EndFor

  // we delete volumes, triangular surfaces and revolution axis
  Delete { Volume{V11,V12,V13,V14}; Surface{TS_1,TS_2,TS_3,TS_4}; Line{A1}; }

  S_1={S_11,S_12,S_13,S_14};

  points~{itp}=P_2;
  itp++;
  points~{itp}=P_3;
  itp++;
  points~{itp}=P_4;
  itp++;
  points~{itp}=P_5;
  itp++;
  points~{itp}=P_6;
  itp++;

  curves~{itc}=E_1;
  If (curves~{itc} < 0)
    curves~{itc}=-E_1;
  EndIf
  itc++;
  curves~{itc}=E_2;
  If (curves~{itc} < 0)
    curves~{itc}=-E_2;
  EndIf
  itc++;
  curves~{itc}=E_3;
  If (curves~{itc} < 0)
    curves~{itc}=-E_3;
  EndIf
  itc++;
  curves~{itc}=E_4;
  If (curves~{itc} < 0)
    curves~{itc}=-E_4;
  EndIf
  itc++;
  curves~{itc}=L_1;
  If (curves~{itc} < 0)
    curves~{itc}=-L_1;
  EndIf
  itc++;
  curves~{itc}=L_2;
  If (curves~{itc} < 0)
    curves~{itc}=-L_2;
  EndIf
  itc++;
  curves~{itc}=L_3;
  If (curves~{itc} < 0)
    curves~{itc}=-L_3;
  EndIf
  itc++;
  curves~{itc}=L_4;
  If (curves~{itc} < 0)
    curves~{itc}=-L_4;
  EndIf
  itc++;

  surfs~{its}=S_1[];
  its++;
  surfs~{its}=S_2;
  its++;
  surfs~{its}=S_3;
  its++;
  surfs~{its}=S_4;
  its++;
  surfs~{its}=S_5;
  its++;

  SL_1=newsl;
  If (LoopsStorage == 1)
    loops~{l}=SL_1;
  EndIf
  Surface Loop(SL_1)={S_1,S_2,S_3,S_4,S_5};
Return
