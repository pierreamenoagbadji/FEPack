// Functions to construct .geo files
General.ExpertMode = 1; // 1 disables the messages for inexperienced users

// Build periodic mesh of rectangle
Function FEPack_Rectangle
  // Points
  P_1 = newp; Point(P_1) = {x1, y1, z1, h0};
  P_2 = newp; Point(P_2) = {x2, y2, z2, h0};
  P_3 = newp; Point(P_3) = {x3, y3, z3, h0};
  P_4 = newp; Point(P_4) = {x4, y4, z4, h0};
  domain_0 = {P_1[], P_2[], P_3[], P_4[]};

  // Lines
  L_1 = newl; Line(L_1) = {P_1, P_2}; Transfinite Line {L_1} = numNodesX;
  L_2 = newl; Line(L_2) = {P_2, P_3}; Transfinite Line {L_2} = numNodesY;
  L_3 = newl; Line(L_3) = {P_3, P_4}; Transfinite Line {L_3} = numNodesX;
  L_4 = newl; Line(L_4) = {P_4, P_1}; Transfinite Line {L_4} = numNodesY;
  domain_1 = {L_1[]};
  domain_2 = {L_2[]};
  domain_3 = {L_3[]};
  domain_4 = {L_4[]};

  // Loops
  LL_1 = newll; Line Loop(LL_1) = {L_1, L_2, L_3, L_4};

  // Plane surfaces
  S_1 = news; Plane Surface(S_1) = {LL_1};

  // Make the mesh structured if specified
  If (is_structured == 1)
    Transfinite Surface {S_1};
  EndIf

  domain_5 = {S_1[]};
Return

// Build periodic mesh of cuboid
Function FEPack_Cuboid
  // Points
  P_1 = newp; Point(P_1) = {x1, y1, z1, h0};
  P_2 = newp; Point(P_2) = {x2, y2, z2, h0};
  P_3 = newp; Point(P_3) = {x3, y3, z3, h0};
  P_4 = newp; Point(P_4) = {x4, y4, z4, h0};
  P_5 = newp; Point(P_5) = {x5, y5, z5, h0};
  P_6 = newp; Point(P_6) = {x6, y6, z6, h0};
  P_7 = newp; Point(P_7) = {x7, y7, z7, h0};
  P_8 = newp; Point(P_8) = {x8, y8, z8, h0};

  // Lines
  L_1  = newl; Line(L_1)  = {P_1, P_2}; Transfinite Line {L_1}  = numNodesX;
  L_2  = newl; Line(L_2)  = {P_2, P_3}; Transfinite Line {L_2}  = numNodesY;
  L_3  = newl; Line(L_3)  = {P_3, P_4}; Transfinite Line {L_3}  = numNodesX;
  L_4  = newl; Line(L_4)  = {P_4, P_1}; Transfinite Line {L_4}  = numNodesY;
  L_5  = newl; Line(L_5)  = {P_5, P_6}; Transfinite Line {L_5}  = numNodesX;
  L_6  = newl; Line(L_6)  = {P_6, P_7}; Transfinite Line {L_6}  = numNodesY;
  L_7  = newl; Line(L_7)  = {P_7, P_8}; Transfinite Line {L_7}  = numNodesX;
  L_8  = newl; Line(L_8)  = {P_8, P_5}; Transfinite Line {L_8}  = numNodesY;
  L_9  = newl; Line(L_9)  = {P_1, P_5}; Transfinite Line {L_9}  = numNodesZ;
  L_10 = newl; Line(L_10) = {P_2, P_6}; Transfinite Line {L_10} = numNodesZ;
  L_11 = newl; Line(L_11) = {P_3, P_7}; Transfinite Line {L_11} = numNodesZ;
  L_12 = newl; Line(L_12) = {P_4, P_8}; Transfinite Line {L_12} = numNodesZ;

  // Loops
  LL_1 = newll; Line Loop(LL_1) = {-L_1, -L_2, -L_3, -L_4};
  LL_2 = newll; Line Loop(LL_2) = {L_5, L_6, L_7, L_8};
  LL_3 = newll; Line Loop(LL_3) = {L_1, L_10, -L_5, -L_9};
  LL_4 = newll; Line Loop(LL_4) = {L_3, L_12, -L_7, -L_11};
  LL_5 = newll; Line Loop(LL_5) = {L_4, L_9, -L_8, -L_12};
  LL_6 = newll; Line Loop(LL_6) = {L_2, L_11, -L_6, -L_10};

  // Plane surfaces with periodic constraints
  S_1 = news; Plane Surface(S_1) = {LL_1};
  S_2 = news; Plane Surface(S_2) = {LL_2};
  Periodic Surface {S_2} = {S_1} Translate {0, 0, z5-z1};

  S_3 = news; Plane Surface(S_3) = {LL_3};
  S_4 = news; Plane Surface(S_4) = {LL_4};
  Periodic Surface {S_4} = {S_3} Translate {0, y3-y1, 0};

  S_5 = news; Plane Surface(S_5) = {LL_5};
  S_6 = news; Plane Surface(S_6) = {LL_6};
  Periodic Surface {S_6} = {S_5} Translate {x2-x1, 0, 0};

  If (is_structured == 1)
    Transfinite Surface {S_1};
    Transfinite Surface {S_2};
    Transfinite Surface {S_3};
    Transfinite Surface {S_4};
    Transfinite Surface {S_5};
    Transfinite Surface {S_6};
  EndIf

  SL_1 = newsl; Surface Loop(SL_1) = {S_1, S_2, S_3, S_4, S_5, S_6};
  domain_1 = {S_6[]};
  domain_2 = {S_5[]};
  domain_3 = {S_4[]};
  domain_4 = {S_3[]};
  domain_5 = {S_2[]};
  domain_6 = {S_1[]};

  // Volume
  V_1 = newv; Volume(V_1) = {SL_1[]};
  If (is_structured == 1)
    Transfinite Volume {V_1};
  EndIf

  domain_7 = {V_1[]};
Return
