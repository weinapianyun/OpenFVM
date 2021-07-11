// Copyright (C) 2004-2005 by OpenFlower Team
// http://sourceforge.net/projects/openflower/
// Licensed under the GNU GPL Version 2.

// Translated from C++ to C by the OpenFVM team in 01/08/2005

#define Back_Bottom_Left   0
#define Back_Bottom_Right  1
#define Front_Bottom_Left  2
#define Front_Bottom_Right 3
#define Back_Top_Left      4
#define Back_Top_Right     5
#define Front_Top_Left     6
#define Front_Top_Right    7

#define EPSILON     1E-6
#define MAXENTITIES 500

typedef struct
{
  double x, y, z;

} oct_data;

typedef struct _oct_node
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double xmid, ymid, zmid;

  int *entities;
  int nbentities;

  int cnodes[8]; 
  struct _oct_node *nodes[8];

} oct_node;

oct_node *root;
oct_node *leafs;

int nbpointers;

int **ipointer;
oct_node **opointer;

int nbleafs;

void OctCreateOctree (double min[3], double max[3], oct_data * Tab, int nbentities);
void OctDestroyOctree ();

