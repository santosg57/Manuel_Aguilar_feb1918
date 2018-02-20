/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "sequence.h"
#include "al.h"
#include "pattern.h"
#include "hit.h"
#include "block.h"
#include "tree.h"
#include "pam_dist.h"


#define PER_WORD 32

static t_node* Node_From_File(FILE** ,t_edge *,t_tree*);




static float Edge_Weight(float w)
/****************************************************************************/
{

	/*
	*/
	return w;


	return w*w;

	return (Pam_Estimated_From_Observed(w));

	/*
	*/
}

t_tree* Tree_From_File(char *filename)
/****************************************************************************/
{
	int c;
	t_node* node;
	t_node* temp_node;
	FILE *file;
	int i;
	t_edge* edge;
	t_tree* tree;

	fprintf(stderr,"Reading tree from file %s\n",filename);

	if ((file=fopen(filename,"r"))==NULL) {
		fprintf(stderr,"file %s could not be opened\n",filename);
		fprintf(stderr,"%s: line %d\n",__FILE__,__LINE__);
		exit (1);
	}

	tree= malloc(sizeof(t_tree));

	node= tree->top_node= (t_node*)calloc(1,sizeof(t_node));
	tree->nr_edges= 0;
	tree->nr_nodes= 0;

	node->index= tree->nr_nodes;
	tree->nodes[ tree->nr_nodes++]= node;

	if ((c= fgetc(file))==EOF)
		return NULL;
	while ((char)c!='(') {
		if ((c= fgetc(file))==EOF)
			return NULL;
	}

	if ((char)c=='(') {
		node->type= intern;
		node->ent.internal= malloc(sizeof(t_internal));

		for (i=0;i<3;i++)
		{
			edge= node->ent.internal->edges[i]= malloc(sizeof(t_edge));
			tree->edges[tree->nr_edges++]= edge;
			edge->L= node;
			edge->R= Node_From_File(&file,edge,tree);
			while (((c=fgetc(file))!=EOF) && (char)c!=':');
			if (c==EOF)
				return NULL;
			if (fscanf(file,"%f",&(edge->dist))<1)
				return NULL;
			edge->dist= Edge_Weight(edge->dist);
			if (i<2) {
				while (((c=fgetc(file))!=EOF) && (char)c!=',');
				if (c==EOF)
					return NULL;
			}
		}
		while (((c=fgetc(file))!=EOF) && (char)c !=')');
	}
	else {
		fprintf (stderr,"Missing first ( in tree file %s\n",filename);
		fprintf (stderr,"%s: line %d\n",__FILE__,__LINE__);
		exit (1);
	}

	return tree;
}

t_node* Node_From_File(FILE** file,t_edge *from_edge,t_tree* tree)
/****************************************************************************/
{
	char c;
	t_node* node;
	t_node* temp_node;
	int i;
	t_edge *edge;

	node= (t_node*)calloc(1,sizeof(t_node));

	node->index= tree->nr_nodes;
	tree->nodes[tree->nr_nodes++]= node;

	c= fgetc(*file);
	while (!isprint(c))
		c= fgetc(*file);

	if (c=='(') {
		node->type= intern;
		node->ent.internal= malloc(sizeof(t_internal));

		for (i=0;i<2;i++) {
			edge= node->ent.internal->edges[i]= malloc(sizeof(t_edge));
			tree->edges[tree->nr_edges++]= edge;
			edge->L= node;
			edge->R= Node_From_File(file,edge,tree);
			while ((c=fgetc(*file))!=':');
			fscanf(*file,"%f",&(edge->dist));
			edge->dist= Edge_Weight(edge->dist);
			if (i<1)
				while ((c=fgetc(*file))!=',');
		}
		while ((c=fgetc(*file))!=')');
		node->ent.internal->edges[2]= from_edge;
	}
	else {
		int len=0;
		node->type= leaf;
		node->ent.leaf= malloc(sizeof(t_leaf));
		node->ent.leaf->name= (char*)calloc(20,sizeof(char));
		node->ent.leaf->name[len++]=c;
		while ((c=fgetc(*file))!=':')
			if (isprint(c))
				node->ent.leaf->name[len++]=c;
		node->ent.leaf->name[len]= '\0';
		node->ent.leaf->edge= from_edge;
		ungetc((int)':',*file);
	}

	return node;
}	


static int Node_Calc_Num(t_node* node,t_node* last)
/****************************************************************************/
{

	if (node->type==leaf) {
		return 1;
	}
	else
	{
		int i, num;
		t_node *n;

		for (i=num=0;i<3;i++){
			n= Edge_Get_Other_Node(node->ent.internal->edges[i],node);
			if (n!=last)
				num+= Node_Calc_Num(n,node);
		}

		return (num);
	}
}

int Tree_Calc_Num(t_tree* tree)
/****************************************************************************/
{
	return Node_Calc_Num(tree->top_node,NULL);
}



static void Node_Print (t_node* node,t_node* last)
/****************************************************************************/
{	
	if (node->type==intern)
	{
		int i;
		t_node *n;
		int started;

		printf ("(");
		started= 0;
		for (i=0;i<3;i++) {
			n= Edge_Get_Other_Node(node->ent.internal->edges[i],node);
			if (n!=last) {
				if (started)
					printf (",");
				Node_Print (n,node);
				printf (":%f",node->ent.internal->edges[i]->dist);
				started++;
			}
		}
		printf (")");
	}
	else
		printf ("%s",node->ent.leaf->name);
}


void Tree_Print (t_tree* tree)
/****************************************************************************/
{	
	printf ("Number of seqs:%d\n",Tree_Calc_Num(tree));
	Node_Print(tree->top_node,NULL);
	printf ("\n");
}


static float Node_Divergence(t_node* node,t_node *last)
/****************************************************************************/
{
	if (node->type==leaf)
		return 0.0;
	else
	{
		int i;
		t_node *n;
		float div;

		div= 0.0;
		for (i=0;i<3;i++) {
			n=Edge_Get_Other_Node(node->ent.internal->edges[i],node);
			if (n!=last) {
				div+= node->ent.internal->edges[i]->dist;
				div+= Node_Divergence(n,node);
			}
		}
		return div;
	}
}

float Tree_Divergence(t_tree *tree)
/****************************************************************************/
{
	float div= 0.0;
	int i;

	for (i=0;i<tree->nr_edges;i++)
		div+= tree->edges[i]->dist;

	printf ("%3d:div:%f\n",tree->nr_edges,div);

	return ( Node_Divergence(tree->top_node,NULL));
}


void Rec_Edge_Set_Of_Nodes(t_node* node, t_node *last, unsigned int *set_of_nodes)
/****************************************************************************/
{
	set_of_nodes[node->index/PER_WORD]|= 1 << (node->index%PER_WORD);

	if (node->type==leaf)
		return;
	else
	{
		int i;
		t_node *n;
		for (i=0;i<3;i++) {
			n=Edge_Get_Other_Node(node->ent.internal->edges[i],node);
			if (n!=last) 
				Rec_Edge_Set_Of_Nodes(n,node,set_of_nodes);
		}
	}
}

void Make_Edge_Sets_Of_Nodes(t_tree *tree)
/****************************************************************************/
{
	int i;
	int j;
	t_edge *edge;
	t_node *node;
	unsigned int *set_of_nodes;
	int index;
	unsigned int *temp_uint;

	tree->node_set_length= tree->nr_nodes/PER_WORD + ((tree->nr_nodes%PER_WORD)>0 ? 1: 0);

	set_of_nodes= calloc(tree->node_set_length,sizeof(unsigned int));

	for (i=0;i<tree->node_set_length-1;i++)
		set_of_nodes[i]= ~0;
	temp_uint= set_of_nodes+tree->node_set_length-1;
	for (i=0;i<tree->nr_nodes%PER_WORD;i++)
		(*temp_uint)|= 1 << (i%PER_WORD);
	/*
		set_of_nodes[tree->node_set_length-1]|= 1 << (i%PER_WORD);
	*/

	for (i=0;i<tree->nr_edges;i++)
	{
		edge= tree->edges[i];
		edge->Left_Nodes= calloc(tree->node_set_length,sizeof(unsigned int));
		Rec_Edge_Set_Of_Nodes(edge->L,edge->R,edge->Left_Nodes);
		edge->Right_Nodes= calloc(tree->node_set_length,sizeof(unsigned int));
		/* 
			Use the property that the left set of nodes is the complement of
			the right set of nodes 
		*/
		for (j=0;j<tree->node_set_length;j++)
			edge->Right_Nodes[j]=set_of_nodes[j]&(~edge->Left_Nodes[j]);
	}

	tree->edge_set_length= tree->nr_edges/PER_WORD + ((tree->nr_edges%PER_WORD)>0 ? 1: 0);

	for (i=0;i<tree->nr_nodes;i++)
	{
		node= tree->nodes[i];		
		index= node->index;

		if (node->type!=leaf)
			continue;

		node->LC= calloc(tree->edge_set_length,sizeof(unsigned int));
		node->RC= calloc(tree->edge_set_length,sizeof(unsigned int));

		assert (node->LC && node->RC);

		for (j=0;j<tree->nr_edges;j++) {
			edge= tree->edges[j];
			if (edge->Left_Nodes[index/PER_WORD]&(1<<(index%PER_WORD)))
				node->LC[j/PER_WORD]|= 1<<(j%PER_WORD);
			if (edge->Right_Nodes[index/PER_WORD]&(1<<(index%PER_WORD)))
				node->RC[j/PER_WORD]|= 1<<(j%PER_WORD);
		}
	}
}


#define MAX_INC 1000

float Find_Divergence_Subset(t_tree *tree, char* filename)
/****************************************************************************/
{
	char *to_be_included[MAX_INC];
	int nr_included;
	int i,j,k;
	FILE *inc_file;
	char buf[100];
	unsigned int *EL, *ER;
	int include;
	t_node *node;
	float div;

	if ((inc_file=fopen(filename,"r"))==NULL) {
		fprintf(stderr,"Find_Divergence_Subset: could not open file %s\n",filename);
		fprintf(stderr,"file %s: line %d\n",__FILE__,__LINE__);
		return -100000.0;
	}

	nr_included= 0;
	while ((fgets(buf,100,inc_file))!=NULL) {
		assert (nr_included<MAX_INC);
		to_be_included[nr_included]= malloc(64);
		sscanf(buf,"%s",to_be_included[nr_included]);
		nr_included++;
	}

	EL= calloc(tree->edge_set_length,sizeof(unsigned int));
	ER= calloc(tree->edge_set_length,sizeof(unsigned int));
	assert (ER && EL);

	for (i=0;i<tree->nr_nodes;i++)
	{
		node= tree->nodes[i];

		if (node->type!=leaf)
			continue;

		for (include=j=0;j<nr_included;j++)
			if (strcmp(node->ent.leaf->name,to_be_included[j])==0)
				include= 1;

		if (include)
		{
			/* include node */
			printf ("   %s\n",node->ent.leaf->name);
			for (j=0;j<tree->edge_set_length;j++) {
				EL[j]|= node->LC[j];
				ER[j]|= node->RC[j];
			}
		}

	}
	/*
	printf ("edges included:\n");
	*/
	div= 0.0;
	for (j=0;j<tree->edge_set_length;j++)
	{
		EL[j]&=ER[j];
		if (EL[j]!=0) {
			for (k=0;k<PER_WORD;k++)
				if (EL[j]&(1<<k))
				{
					printf ("   %f\n",tree->edges[j*PER_WORD+k]->dist);
					div+= tree->edges[j*PER_WORD+k]->dist;
				}
		}
	}

	free (ER);
	free (EL);

/*
	printf ("Resulting divergence:%f\n",div);
*/
	printf ("\n");

	return div;
}

int Tree_Make_Seq_Pointers(t_tree *tree, t_sequence **seqs, int nr)
/****************************************************************************/
{
	int i,j;
	int index;
	char *name;

	for (i=0;i<nr;i++) 
		seqs[i]->tree_index= (-1);
		
	for (i=0;i<tree->nr_nodes;i++) {

		if (tree->nodes[i]->type!=leaf)
			continue;

		name= tree->nodes[i]->ent.leaf->name;

		index= -1;
		for (j=0;j<nr;j++)
			if (strcmp(name,seqs[j]->filename)==0)
				index= j;

		if (index<0)
			return 0;
		else {
			tree->nodes[i]->ent.leaf->seq_index= index;
			seqs[index]->tree_index= i;
		}
	}

	for (i=0;i<nr;i++) 
		if (seqs[i]->tree_index==(-1))
			return 0;

	return 1;
}

float Divergence(t_tree *tree,int Inc[] ,t_sequence **seqs,int nr)
/****************************************************************************/
{
	unsigned int *leaf_included;
	unsigned int *EL, *ER;
	int i,j,k;
	int index;
	float div;
	t_node *node;

	leaf_included=calloc(1+tree->nr_nodes/PER_WORD,sizeof(unsigned int));

	for (i=0;i<nr;i++) {
	/*
		j= block->start[i];
		while (hits[j]==0 && j<block->start[i+1])
			j++;
		if (j<block->start[i+1]) 
	*/

		if (Inc[i]==1){
			index= seqs[i]->tree_index;
			leaf_included[index/PER_WORD]|=1<<(index%PER_WORD);
		}
	}

	EL= calloc(tree->edge_set_length,sizeof(unsigned int));
	ER= calloc(tree->edge_set_length,sizeof(unsigned int));
	assert (ER && EL);

	for (i=0;i<tree->nr_nodes;i++)
	{
		node= tree->nodes[i];

		if (node->type!=leaf)
			continue;

		if (leaf_included[i/PER_WORD]&(1<<i%PER_WORD))
		{
			/* include node */
			/*
			printf ("   %s\n",node->ent.leaf->name);
			*/
			for (j=0;j<tree->edge_set_length;j++) {
				EL[j]|= node->LC[j];
				ER[j]|= node->RC[j];
			}
		}

	}
	/*
	printf ("edges included:\n");
	*/
	div= 0.0;
	for (j=0;j<tree->edge_set_length;j++)
	{
		EL[j]&=ER[j];
		if (EL[j]!=0) {
			for (k=0;k<PER_WORD;k++)
				if (EL[j]&(1<<k))
				{
				/*
					printf ("   %f\n",tree->edges[j*PER_WORD+k]->dist);
				*/
					div+= tree->edges[j*PER_WORD+k]->dist;
				}
		}
	}

	free (ER);
	free (EL);

/*
	printf ("Resulting divergence:%f\n",div);
*/

	free (leaf_included);

	return div;
}
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
