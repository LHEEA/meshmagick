h = DefineNumber[ 10, Name "Parameters/hauteur" ];
r = DefineNumber[ 2.5, Name "Parameters/rayon" ];
elt_size = DefineNumber[ 0.1, Name "Parameters/element size" ];


Point(1) = {0, 0, -h/2, elt_size};
Point(2) = {r, 0, -h/2, elt_size};
Point(3) = {0, r, -h/2, elt_size};

Circle(1) = {3, 1, 2};
Line(2) = {2, 1};
Line(3) = {1, 3};

Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};

Extrude {0, 0, h/2} {
  Line{1};
}

Transfinite Surface {5} Alternated;
Recombine Surface {5}; // A decommenter si tu veux des quadrangles aux extremites

Transfinite Surface {9} Alternated;
Recombine Surface {9}; // A decommenter si tu veux des quadrangles autour
