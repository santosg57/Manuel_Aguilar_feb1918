/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef TREE_DEF
#define TREE_DEF


typedef enum {intern,leaf} t_node_type;

typedef struct 
{
	struct t_edge_temp *edges[3]; 
} t_internal;

typedef struct
{
	struct t_edge_temp *edge;
	char* name;
	int   index;
	int   order;
	int   seq_index;
} t_leaf;

typedef struct t_node_temp
{
	int index;
	unsigned int *LC;
	unsigned int *RC;
	t_node_type type;
	union  { t_internal* internal; t_leaf* leaf;} ent;
} t_node;

typedef struct t_edge_temp
{
	float dist;
	t_node *L, *R;
	unsigned int *Left_Nodes;
	unsigned int *Right_Nodes;
} t_edge;

#define Edge_Get_Other_Node(e,n) (((e)->L==(n))?(e)->R:(e)->L)

typedef struct 
{
	t_edge *edges[1000];
	t_node *nodes[1000];
	int nr_edges;
	int nr_nodes;
	int node_set_length;
	int edge_set_length;
	t_node *top_node;
} t_tree;
	 


t_tree* Tree_From_File(char *);
void Tree_Print (t_tree*);
float Tree_Divergence(t_tree*);
void Make_Edge_Sets_Of_Nodes(t_tree *);
void Make_Random_Sets_With_Divergence(t_tree *, float );
float Find_Divergence_Subset(t_tree*,char*);

int Tree_Make_Seq_Pointers(t_tree*, t_sequence**, int );
float Divergence(t_tree*,int[],t_sequence**,int);

#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
