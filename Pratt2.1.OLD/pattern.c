/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "menu.h"
#include "hit.h"


#include "sequence.h"
#include "al.h"
#define PAT_MAIN
#include "pattern.h"
#include "block.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

/*
#define DEBUG
*/
extern unsigned int Member_In_Set[];
extern int Number_In_Set[];
extern float Freq_Set[];
extern float Bit_Set[];
extern int Nr_Sets;
extern int Nr_Sets_Block;
extern int Nr_Sets_Unit;
extern int Nr_Sets_In_Search;
extern float Shannon;
extern int Num_Ones[];

extern  t_Block_Options *B_Options;
extern t_Generalize_Symbol_List **Gen_Symbol_Lists;

extern int Length_Block;
extern int nrseqs;

#ifdef DEBUG
#define CHECK(a) assert(a)
#else
#define CHECK(a)
#endif

void Pat_Print_gcg(FILE* file,t_Pat_Info *pat)
/******************************************************************/
{
	int i,j,k;
	int symbol;
	unsigned int mask;
	int num;

	for (i=0;i<pat->num_comp;i++) 
	{
		if (pat->gaps[i]+pat->flexi[i]>0)
			fprintf (file,"x{%d,%d}",pat->gaps[i],pat->gaps[i]+pat->flexi[i]);
		symbol= pat->positions[i];
		if (symbol>=0)
		{
			if (Number_In_Set[symbol]>1)
				fprintf (file,"("); 
			num= 0;
			for (j=0;j<26;j++)
				if (Member_In_Set[symbol]&(1<<j))
				{
					if (num>0)
						fprintf(file,",");
					fprintf (file,"%c",(char)(j+'A'));
					num++;
				}
			if (Number_In_Set[symbol]>1) 
				fprintf (file,")"); 
		}
		else
		{
			assert (symbol==(-1));
			fprintf (file,"(");
			mask= (unsigned int) (pat->spec_pos[i]);
			num= 0;
			for (j=0;j<Nr_Sets_Unit;j++) {
				if ((mask)&(1<<j)) {
					for (k=0;k<26;k++)
						if (Member_In_Set[j]&(1<<k))
						{
							if (num>0)
								fprintf (file,",");
							fprintf (file,"%c",(char)(k+'A'));
							num++;
						}
				}
			}
			fprintf (file,")");
		}
	}
/*
	fprintf (file,"\n");
*/
}


static void Pattern_Print(FILE* file,t_Pat_Info *pat)
/******************************************************************/
{
	int i,j;
	int more;
	unsigned int pos;

	for (i=0;i<pat->num_comp;i++) 
	{
		for (j=0;j<pat->gaps[i];j++)
			fprintf (file,"x");
		for (j=0;j<pat->flexi[i];j++)
			fprintf (file,"-");

		pos= pat->spec_pos[i];

		if ((more=NUM_ONES32(pos))>1)
			fprintf (file,"["); 
		for (j=0;j<26;j++)
			if (pos&(1<<j))
				fprintf (file,"%c",(char)(j+'A'));
		if (more>1)
			fprintf (file,"]"); 
	}
	fprintf (file,"\n");
}

void Pat_Print(FILE* file,t_Pat_Info *pat)
/******************************************************************/
{
	int i,j;
	int more;
	unsigned int pos;

	if (!B_Options->Prosite_Style) {
		Pattern_Print (file,pat);
		return;
	}

	for (i=0;i<pat->num_comp;i++) 
	{
		pos= pat->spec_pos[i];
		if (i>0)
			fprintf (file,"-");
		if ((pat->gaps[i]+pat->nflexi[i])>0)
		{
			fprintf (file,"x");
			if (pat->gaps[i]!=1 || pat->nflexi[i]>0)
			{	
				if (pat->nflexi[i]>0)
					fprintf(file,"(%d,%d)",pat->gaps[i],pat->gaps[i]+pat->nflexi[i]);
				else
					fprintf(file,"(%d)",pat->gaps[i]);
			}
			fprintf (file,"-");
		}

		if ((more=NUM_ONES32(pos))>1)
			fprintf (file,"["); 
		for (j=0;j<26;j++)
			if (pos&(1<<j))
				fprintf (file,"%c",(char)(j+'A'));
		if (more>1)
			fprintf (file,"]"); 
	}
/*
	fprintf (file,"\n");
*/
}

/*
void Pat_Print_Site(FILE* file,t_Pat_Info *pat)
{
	int i,j,k;
	int symbol;
	unsigned int mask;

	for (i=0;i<pat->num_comp;i++) 
	{
		symbol= pat->positions[i];
		if (i>0)
			fprintf (file,"-");
		if ((pat->gaps[i]+pat->nflexi[i])>0)
		{
			fprintf (file,"x");
			if (pat->gaps[i]!=1 || pat->nflexi[i]>0)
			{	
				if (pat->nflexi[i]>0)
					fprintf(file,"(%d,%d)",pat->gaps[i],pat->gaps[i]+pat->nflexi[i]);
				else
					fprintf(file,"(%d)",pat->gaps[i]);
			}
			fprintf (file,"-");
		}

		if (symbol>=0)
		{
			fprintf (file,"%c",(char)symbol);
		}
		else
		{
			assert (symbol==(-1));
			fprintf (file,"[");
			mask= (unsigned int) (pat->spec_pos[i]);
			for (j=0;j<26;j++) {
				if ((mask)&(1<<j)) {
					fprintf (file,"%c",(char)(j+'A'));
				}
			}
			fprintf (file,"]");
		}
	}
}
*/

int Pat_Copy (t_Pattern* p1, t_Pattern* p2)
/******************************************************************/
{
	int i;
	int n;
	size_t length;

	CHECK (p1 && p1->debug==PDEBUG);
	CHECK (p2 && p2->debug==PDEBUG);

	p2->length= p1->length;
	p2->num_comp= p1->num_comp;
	p2->info= p1->info;
	p2->num_flex= p1->num_flex;

	CHECK (p2->positions && p2->gaps);

	n= p1->num_comp;
	length= n*sizeof(int);
	memcpy(p2->positions,p1->positions,length);
	memcpy(p2->gaps,p1->gaps,length);
	memcpy(p2->flexi,p1->flexi,length);

	memcpy(p2->spec_pos,p1->spec_pos,n*sizeof(unsigned int));

	p2->num_sub= p1->num_sub;
	p2->num_sub_new= p1->num_sub_new;

	memcpy (p2->lengths,p1->lengths,p1->num_sub*sizeof(int));

	CHECK (p1->debug==PDEBUG);
	CHECK (p2->debug==PDEBUG);

	return 1;
}



void Pat_Check(t_Pattern *p1)
/******************************************************************/
{
	int i;

	assert (p1->debug== PDEBUG);
	assert (p1->num_comp<=B_Options->Max_Num_Comp);
	assert (p1->positions);
	assert (p1->gaps);
	assert (p1->num_sub<=B_Options->Max_Flex_Prod);
	assert (p1->lengths);
	assert (p1->bit_vecs);
	assert (p1->flexi);
	for (i=0;i<p1->num_sub;i++)
		assert (p1->bit_vecs[i]);
}	

	

t_Pat_Info* Pat_Info_From_Pattern(t_Pattern* Pattern, int flag)
/******************************************************************/
{
	t_Pat_Info *pat;
	int j,k;
	unsigned int* vec1;
	unsigned int* vec2;
	unsigned int* end2;
	unsigned int* pt2;
	size_t num_bytes;

	pat= malloc(sizeof(t_Pat_Info));
	assert (pat);
	pat->num_comp= Pattern->num_comp;
	pat->length= Pattern->length;
	pat->positions= malloc(B_Options->Max_Num_Comp*sizeof(int));
	pat->spec_pos= malloc(B_Options->Max_Num_Comp*sizeof(unsigned int));
	pat->gaps= malloc(B_Options->Max_Num_Comp*sizeof(int));
	pat->flexi= malloc(B_Options->Max_Num_Comp*sizeof(int));
	pat->nflexi= malloc(B_Options->Max_Num_Comp*sizeof(int));
	pat->used= calloc(B_Options->Max_Num_Flex,sizeof(unsigned int));
	assert (pat->positions && pat->spec_pos && pat->gaps && pat->flexi && pat->nflexi );

	num_bytes= Pattern->num_comp*sizeof(int);
	memcpy(pat->positions,Pattern->positions,num_bytes);
	memcpy(pat->gaps,Pattern->gaps,num_bytes);
	memcpy(pat->flexi,Pattern->flexi,num_bytes);
	memcpy(pat->nflexi,Pattern->flexi,num_bytes);
	memcpy(pat->spec_pos,Pattern->spec_pos,Pattern->num_comp*sizeof(unsigned int));

	num_bytes= Length_Block*sizeof(unsigned int);
	vec2=pat->bit_vector= calloc(Length_Block,sizeof(unsigned int));
	assert (vec2);
	for (j=0; j<Pattern->num_sub;j++) {
		vec1= Pattern->bit_vecs[j];
		pt2= vec2;
		end2= vec2+Length_Block;
		while (pt2<end2)
			*pt2++ |= *vec1++;
	}

	pat->num_seq= Pattern->num_seq;
	pat->information= Pattern->info;
	pat->info_MDL= Pattern->info_MDL;
	pat->nr_match_swiss= pat->nr_occ_swiss= 0;

	pat->special_flag= 0;

	if (flag) {
		if (B_Options->Diagnostic) {
			pat->matching_strings= malloc(MAX_MATCHING_STRINGS*sizeof(char*));
			pat->seq_indexes= malloc(MAX_MATCHING_STRINGS*sizeof(int));
			pat->nr_matching_strings= 0;
			pat->max_nr_matching_strings= MAX_MATCHING_STRINGS;
			assert (pat->matching_strings);
		}

/*
		pat->num_sub= Pattern->num_sub;
		pat->bit_vecs= malloc(B_Options->Max_Flex_Prod*sizeof(unsigned int*));
		assert (pat->bit_vecs);
		for (j=0;j<Pattern->num_sub;j++) {
			vec2= pat->bit_vecs[j]= malloc(num_bytes);
			assert (vec2);
			memcpy(vec2,Pattern->bit_vecs[j],num_bytes);
		}
*/
	}
	else { 	
	/*
		pat->num_sub= 0;
		pat->bit_vecs= NULL;
	*/
	}

	/* 
		Inge Sept 10th -- moved this out of if (flag) in order to 
		store the bitvectors also after refinement  -- these are
		needed to print out the patterns in new format
	*/
	pat->num_sub= Pattern->num_sub;
	pat->bit_vecs= malloc(B_Options->Max_Flex_Prod*sizeof(unsigned int*));
	assert (pat->bit_vecs);
	for (j=0;j<Pattern->num_sub;j++) {
		vec2= pat->bit_vecs[j]= malloc(num_bytes);
		assert (vec2);
		memcpy(vec2,Pattern->bit_vecs[j],num_bytes);
	}
	/* end of what I moved out of if (flag) */

	return pat;
}


int Pat_Info_From_String(char *string, t_Pat_Info *pat)
/******************************************************************/
{
	int i,lower,upper;

	pat->length= 0;
	pat->num_comp= 0;
	pat->gaps[0]= 0;
	pat->flexi[0]= 0;
	pat->num_comp= 0;

	i= 0;
	while (i<(int)strlen(string)-1)
	{
		if (string[i]=='x')
		{
			if (string[i+1]=='-') {
				pat->gaps[pat->num_comp]= 1;
				pat->flexi[pat->num_comp]= 0;
			}
			else if (strstr(string+i+1,",")!=NULL && 
						strstr(string+i+1,",")<strstr(string+i+1,")")) {
				sscanf(string+i+1,"(%d,%d)",&lower,&upper); 
				pat->gaps[pat->num_comp]= lower;
				pat->flexi[pat->num_comp]= upper-lower;
			}
			else if (sscanf(string+i+1,"(%d)",&lower)) {
				pat->gaps[pat->num_comp]= lower;
				pat->flexi[pat->num_comp]= 0;
			}
			else {
				printf ("wildcard regions non-parsable\n");
				return 0;
			}
			while (string[i]!='-' && string[i]!='\n') 
				i++;
			i++;
		}
		else
		{
			pat->gaps[pat->num_comp]= pat->flexi[pat->num_comp]= 0;
		}

		if (isupper(string[i])) 
		{
			pat->positions[pat->num_comp]=string[i];
			i++;
			assert (string[i]=='-' || string[i]=='\n');
			i++;
		}
		else if (string[i]=='[')
		{
			i++;
			pat->positions[pat->num_comp]= (-1);
			pat->spec_pos[pat->num_comp]= 0;
			while (string[i]!=']') {
				if (isupper(string[i]))
					pat->spec_pos[pat->num_comp]|= 1<<(string[i]-'A');
				else {
					printf ("error:%c not upper case (within brackets)\n",string[i]);
					return 0;
				}
				i++;
			}
			i++;
			assert (string[i]=='-' || string[i]=='\n');
			i++;
		}
		pat->length+= 1+pat->gaps[pat->num_comp]+pat->flexi[pat->num_comp];
		pat->num_comp++;
	}
	return 1;
}


void Pat_Info_Free(t_Pat_Info* pat)
/******************************************************************/
{
	int i;

	return;

	if (pat==NULL)
		return;
	if (pat->num_comp<0 || pat->num_comp>1000)
		return;

	if (pat->positions)
		free (pat->positions);
	if (pat->spec_pos)
		free (pat->spec_pos);
	if (pat->gaps)
		free (pat->gaps);
	if (pat->flexi)
		free (pat->flexi);
	if (pat->bit_vector)
		free (pat->flexi);
	if (pat->bit_vecs) {
		for (i=0;i<pat->num_sub;i++)
			if (pat->bit_vecs[i])
				free (pat->bit_vecs[i]);	
		free (pat->bit_vecs);
	}
	free (pat);
}


int Pattern_Equal_Pat_Info(t_Pattern* Pat, t_Pat_Info *Pat_Info)
/******************************************************************/
{
	int i;

	if (Pat->num_comp!=Pat_Info->num_comp) 
		return 0;

	for (i=0;i<Pat->num_comp;i++) {
		if ((Pat->positions[i]!=Pat_Info->positions[i]) ||
			 (Pat->gaps[i]!=Pat_Info->gaps[i]) ||
			 (Pat->flexi[i]!=Pat_Info->flexi[i]) )
			return 0;
	}
	return 1;
}

static int Pat_Info_Same_Pattern(t_Pat_Info* P1, t_Pat_Info* P2)
/******************************************************************/
{
	int i;
	int same;

	same= 1;

	if (P1->num_comp!=P2->num_comp)
		same= 0;

	for (i=0;(same && i<P1->num_comp);i++)
	{
		if (P1->gaps[i]!=P2->gaps[i] || P1->nflexi[i]!=P2->nflexi[i])
			same= 0;

		if (P1->spec_pos[i]!=P2->spec_pos[i])
			same= 0;
	}

	return same;
}

int Pat_Info_Match_Same_Segments(t_Pat_Info* P1, t_Pat_Info* P2)
/******************************************************************/
/* check if the set of segments matched by P2 is also matched by  */
/* P1. P1 scores higher than P2, and if they match the same       */
/* segments, P2 is not reported                                   */
/* 21/6/95 I. Jonassen                                            */
/******************************************************************/
{
	int i,j,k;
	int l2, l1,diff;
	int extra;
	int found;
	unsigned int *b1, *b2;
	int pos;

	if (P1==NULL || P2==NULL)
		return 0;

	if (Pat_Info_Same_Pattern(P1,P2))
		return 1;

	b1= P1->bit_vector;
	b2= P2->bit_vector;
	l1= P1->length; 
	l2= P2->length;
	diff= l1-l2;

	if (diff<0)
		return 0;

	extra= B_Options->Max_Flex * B_Options->Max_Num_Flex;

	if (extra>0)
		return 0;

	for (i=0;i<Length_Block;i++)
	{
		if (b2[i]!=0) 
		{
			for (j=0;j<PER_WORD;j++)
			{
				if (b2[i]&(1<<j))
				{
					pos= i*PER_WORD+j;
					for (found=0,k=pos-diff-extra;k<=pos+extra;k++) {
						if (b1[k/PER_WORD]&(1<<(k%PER_WORD)))
							found= 1;
					}
					if (found==0)
						return 0;
				}
			}
		}
	}
	return 1;
}

int Pat_Info_Num_X(t_Pat_Info *pat)
/******************************************************************/
{
	int i, num;

	for (i=num=0;i<pat->num_comp;i++)
	{
		if (pat->gaps[i]+pat->flexi[i]>0)
			num++;
	}
	return num;
}


int Pat_Info_Num_Char(t_Pat_Info *pat)
/******************************************************************/
{
	int i, num;

	for (i=num=0;i<pat->num_comp;i++)
	{
		num++;
		/*
		if (pat->gaps[i]+pat->flexi[i]>0)
			num++;
		*/
	}
	return num;
}

int Pattern_Num_X(t_Pattern *pat)
/******************************************************************/
{
	int i, num;

	for (i=num=0;i<pat->num_comp;i++)
	{
		if (pat->gaps[i]+pat->flexi[i]>0)
			num++;
	}
	return num;
}


int Pattern_Num_Char(t_Pattern *pat)
/******************************************************************/
{
	int i, num;

	for (i=num=0;i<pat->num_comp;i++)
	{
		num++;
		/*
		if (pat->gaps[i]+pat->flexi[i]>0)
			num++;
		*/
	}
	return num;
}



t_Pat_Lead_Info* Lead_Info_From_Pat(t_Pattern *pat,t_Pat_Lead_Info *lead)
/******************************************************************/
{
	t_Pat_Lead_Info *best;
	int j;

	if (lead)
		best= lead;
	else {
      best= malloc(sizeof(t_Pat_Lead_Info));
      assert (best);
      best->positions= malloc(B_Options->Max_Num_Comp*sizeof(int));
      best->gaps= malloc(B_Options->Max_Num_Comp*sizeof(int));
      best->flexi= malloc(B_Options->Max_Num_Comp*sizeof(int));
      assert (best->positions && best->gaps && best->flexi);
   }

   best->score= pat->info;
   best->num_comp= pat->num_comp;
   best->num_sub= pat->num_sub;
   best->num_flex= pat->num_flex;
   for (j=0;j<pat->num_comp;j++) {
      best->positions[j]= pat->positions[j];
      best->gaps[j]= pat->gaps[j];
      best->flexi[j]= pat->flexi[j];
   }

	return best;
}



/******************************************************************/
/********** ROUTINES ONLY USED LOCALLY                   **********/
/******************************************************************/
static void Pat_Graph_Make_Edges(t_Pat_Graph *);
static void Pat_Graph_Find_Length_Paths(t_Pat_Graph *);
static void Graph_Calc_Upper_Bounds(t_Pat_Graph *);
static void Possible_Edges_In_Graph(t_Pat_Graph*,t_block*);
/******************************************************************/

t_Pat_Graph* Pat_Graph_From_Alignment(t_alignment* al,t_block *block)
/******************************************************************/
{
	t_Pat_Graph *temp;
	int i,j,k;
	unsigned int column;
	int *distance;
	int hopeless;
	int symbol;
	int *position;
	int last;

	printf ("Constructing and analysing pattern graph that\n");
	printf ("  is used to guide the initial pattern search\n");

	temp= malloc(sizeof(t_Pat_Graph));
	assert (temp);

	assert (B_Options->Max_Num_Comp<G_MAX_COMP && 
			  B_Options->Max_Num_Flex<G_MAX_FLEX);

	temp->alignment= al;
	temp->nr_seq_alignment= al->nr_seq;
	temp->max_nr_nodes= al->nr_pos;
	temp->nr_nodes= 0;
	assert (temp->max_nr_nodes>0);
	temp->nodes= malloc(temp->max_nr_nodes*sizeof(t_Pat_Node*));
	assert (temp->nodes);

	position= calloc(al->nr_seq,sizeof(int));
	/*
	printf ("al->nr_pos=%d\n",al->nr_pos);
	*/
	last= 0;

	for (i=0;i<al->nr_pos;i++)
	{
		column= 0;
		hopeless= 0;
		for (j=0;j<al->nr_seq;j++) {
			if (isupper(al->positions[j][i])) {
				position[j]++;
				column|= 1<<(al->positions[j][i]-'A');
				if ((al->seq_index) && (al->seq_index[i]!=last)) {
					position[j]+= B_Options->Max_Gap+1;
					last= al->seq_index[i];
				}
			}
			else {
				hopeless= 1;
			}
		}
		if (!hopeless)
		{
			symbol= (-1); j= 0;
			while (symbol<0 && j<Nr_Sets) {
				if ((column & Member_In_Set[j])==column)
					symbol= j;
				j++;
			}
		}
		else
			symbol= (-1);

/*
		if (symbol>=0) {
*/
		if (symbol>=0 && symbol<B_Options->Nr_Symbol_Search1) {
			t_Pat_Node *Node;

			Node= temp->nodes[temp->nr_nodes]= malloc(sizeof(t_Pat_Node));	
			assert (Node);

			Node->seq_index= (al->seq_index==NULL ? 0 : al->seq_index[i]);

			Node->symbol= symbol;
			Node->edges= malloc((B_Options->Max_Gap+1)*sizeof(t_Pat_Edge*));
			Node->max_nr_edges= B_Options->Max_Gap+1;
			Node->nr_edges= 0;
			Node->pos_alignment=malloc(sizeof(int)*temp->nr_seq_alignment);
			assert (Node->pos_alignment);

			for (j=0;j<temp->nr_seq_alignment;j++)
				Node->pos_alignment[j]= position[j];
			temp->nr_nodes++;
		}
	}

/*
	printf ("temp->nr_nodes=%d\n",temp->nr_nodes);
*/
	free (position);

	Pat_Graph_Make_Edges(temp);
	Pat_Graph_Find_Length_Paths(temp);
	Possible_Edges_In_Graph(temp,block);
	Graph_Calc_Upper_Bounds(temp);

	printf ("-- finished constructing pattern graph.\n\n");
	
	return temp;
}


static void Pat_Graph_Make_Edges(t_Pat_Graph *graph)
/******************************************************************/
{
	int i,j,k;
	int max_diff_allowed;
	int max_flex_allowed;
	int min_gap,max_gap;
	int gap;
	t_Pat_Edge *edge;
	t_Pat_Node *node1, *node2;

	max_diff_allowed= B_Options->Max_Gap+1;
	max_flex_allowed= B_Options->Max_Flex;

	for (i=0;i<graph->nr_nodes;i++)
	{
		node1= graph->nodes[i];
		j= i+1;
		while ((j<graph->nr_nodes) && 
				 (graph->nodes[j]->pos_alignment[0]<=node1->pos_alignment[0]+max_diff_allowed))
		{
			node2= graph->nodes[j];
			min_gap= 1000; max_gap= 0;
			for (k=0;k<graph->nr_seq_alignment;k++) {
				gap= node2->pos_alignment[k]-node1->pos_alignment[k];
				min_gap=min(min_gap,gap);
				max_gap=max(max_gap,gap);
			}
			if ((max_gap<=max_diff_allowed) && ((max_gap-min_gap)<=max_flex_allowed))
			{
				edge= node1->edges[node1->nr_edges++]= malloc(sizeof(t_Pat_Edge));	
				assert (edge);
				edge->min_aa= min_gap-1;
				edge->max_aa= max_gap-1;
				edge->Node= node2;
			}
			j++;
		}
	}
}


static void Graph_Calc_Upper_Bounds(t_Pat_Graph *graph)
/******************************************************************/
{
	int i,j,k;
	t_Pat_Node *node,*node2;
	t_Pat_Edge *edge;
	int nr_edges;
	int c,n;
	float info;
	float red;
	float *fpt,*ept;
	int f,flex;
	
	for (i=graph->nr_nodes-1;i>=0;i--) {
		node= graph->nodes[i];
		info= Bit_Set[node->symbol];
		fpt=&(node->Max_Score[0][0]);
		ept=&(node->Max_Score[G_MAX_COMP-1][G_MAX_FLEX-1]);
		while (fpt<=ept)
			*fpt++= info;
		for (j=0;j<node->nr_edges;j++) {
			edge=node->edges[j];
			flex= edge->max_aa - edge->min_aa;
			if (edge->possible[0][0]==0)
				continue;
			else if (((edge->possible[0][0]&(1<<16))==0) && (flex==0))
				flex++;
			red= (float)flex*0.5;
			f= (flex > 0 ? 1 : 0);
			node2= edge->Node;
			for (c=1;c<G_MAX_COMP;c++)
				for (n=f;n<G_MAX_FLEX;n++)
					node->Max_Score[c][n]= max(node->Max_Score[c][n],
				   	                    node2->Max_Score[c-1][n-f]+info-red);
		}
		/*
		for (c=0;c<G_MAX_COMP;c++) {
			for (n=0;n<G_MAX_FLEX;n++)
				printf ("%5.1f ",node->Max_Score[c][n]);
			printf ("\n");
		}
		*/
	}
}


static void Pat_Graph_Find_Length_Paths(t_Pat_Graph *graph)
/******************************************************************/
{
	int i,j,k;
	t_Pat_Node *node1, *node2;
	t_Pat_Edge *edge;
	int extension;
	float max_extra_score;

	for (i=graph->nr_nodes-1;i>=0;i--) {
		node1= graph->nodes[i];
		extension= 0;
		max_extra_score= 0.0;
		for (j=0;j<node1->nr_edges;j++) {
			edge= node1->edges[j];
			node2= edge->Node;
			extension= max(extension,node2->length_longest_path);
			max_extra_score= max(max_extra_score,
					 (node2->max_score-0.5*(float)(edge->max_aa-edge->min_aa)));
		}
		node1->length_longest_path= 1+extension;
		node1->max_score= Bit_Set[node1->symbol]+max_extra_score;

		/*
		printf ("length=%3d, score= %8.4f\n",node1->length_longest_path,node1->max_score);
		*/
	}
}



static void Possible_Edges_In_Graph(t_Pat_Graph *graph, t_block *block)
/******************************************************************/
{
	int i,j,k,l;
	t_Pat_Node *node1, *node2;
	t_Pat_Edge *edge;
	int symbol1,symbol2;
	int i1,i2;
	int s1,s2;
	t_Generalize_Symbol_List *l1, *l2;
	int flex;
	int possible;
	int max_flex;

	max_flex=B_Options->Max_Flex;

	/*
	printf ("Possible_Edges_In_Graph\n");
	*/
	Matching_Prepare(nrseqs,max_flex);

	for (i=0;i<graph->nr_nodes;i++) {
		node1= graph->nodes[i];
		for (j=0;j<node1->nr_edges;j++) {
			edge= node1->edges[j];
			flex= edge->max_aa-edge->min_aa;
			node2= edge->Node;
			symbol1= node1->symbol; l1= Gen_Symbol_Lists[symbol1];
			symbol2= node2->symbol; l2= Gen_Symbol_Lists[symbol2];

			for (i1=0;i1<l1->nr_symbols;i1++)
				for (i2=0;i2<l1->nr_symbols;i2++)
					edge->possible[i1][i2]= 0;

			for (i1=0;i1<l1->nr_symbols;i1++)
			{
				s1= l1->symbols[i1];
				for (i2=0;i2<l2->nr_symbols;i2++)
				{
					s2= l2->symbols[i2];
					edge->possible[i1][i2]= Nr_Match(s1,s2,edge->min_aa,edge->max_aa,max_flex,block);
				}
			}
		}
	}
	/*
	printf ("end Possible_Edges_In_Graph\n");
	*/
}

static int Match_To_String (t_Pat_Info *Pat_Info, int i1, int vector_index, char string[],t_sequence* seq,int flag)
/******************************************************************/
{
	int i,j;
	int len;
	int l;
	int index;
	int prod;
	int extra;
	int fi; /* index to flexibility */

	index=i1; 
	prod= 1;
	len= l= 0;

	fi= 0;
	for (i=0;i<Pat_Info->num_comp;i++) {
		for (j=0;j<Pat_Info->gaps[i];j++) {
			string[len++]= tolower(seq->seq[index++]);
			l++;
		}
		if (Pat_Info->flexi[i]>0) {
		/*
			for (j=0;j<((vector_index/prod)%(prod*(1+Pat_Info->flexi[i])));j++) {
		*/
			extra= ((vector_index/prod)%(1+Pat_Info->flexi[i]));
			Pat_Info->used[fi]|= 1<<extra;
			for (j=0;j<extra;j++) {
				string[len++]= tolower(seq->seq[index++]);
				l++;
			}
			for (    ;j<Pat_Info->nflexi[i];j++) {
				string[len++]= '-';
			}
			prod*=(Pat_Info->flexi[i]+1); 
			fi++;
		}
		string[len++]= toupper(seq->seq[index++]);
		l++;
	}
	string[len]='\0';

	return l;
}

void Special_Print_Pattern_Matches_For_One_Sequence(FILE* file,void* block_pointer,t_sequence** seqs, int nr, t_Pat_Info* Pat_Info)
/******************************************************************/
{
   unsigned int* stat;
   int length;
   int i,j,k,l;
	unsigned int filter;
	t_sequence *seq;
	int i1,i2;
	int ii;
	t_block* block;
	char match_string[1000];
	int len;

	block= (t_block*)block_pointer;
 
   stat= Pat_Info->bit_vector;
   length= Pat_Info->length;
 
	if (file) 
		fprintf (file,"Occurrences: %d(%d)\n",Pat_Info->nr_occ,Pat_Info->nr_match);
 
	seq= seqs[0];

   for (j=0;j<block->len;j++) {
      for (k=0;k<PER_WORD;k++) {
			filter= (1<<k);
         if (stat[j] & filter ) {
				for (l=0;l<Pat_Info->num_sub;l++) {
					if (Pat_Info->bit_vecs[l][j] & filter) {

						i1= (int)(j*PER_WORD+k); 
						len= Match_To_String (Pat_Info,i1,l,match_string,seq,(file!=NULL));
						i2= i1+len-1;
						if (file) {
   						fprintf (file,"%14s%c:%6d-%6d: ",seq->filename,(seq->restricted?'*':' '),i1+1,i2+1);
							for (ii= i1-5;ii<i1;ii++)
      				 		fprintf (file,"%c",(ii>=0 ? (tolower(seq_element(seq,ii))) : ' ') );

							fprintf (file," %s ",match_string);
   						for (ii=i2+1;ii<=i2+5;ii++)
      						fprintf (file,"%c",(ii<seq->length ? (tolower(seq_element(seq,ii))) : ' ') );
   						fprintf(file,"\n");
						}
					}
				}
			}
      }
	}
}


void Special_Print_Pattern_Matches(FILE* file,void* block_pointer,t_sequence** seqs, int nr, t_Pat_Info* Pat_Info)
/******************************************************************/
/*
	Added by Inge Sept 10th 1996 to get pattern matches printed out aligned
	(with gaps) and with capital letters matching pattern symbols and lower
	case letters matching wildcards

	Modified Jan 3rd 1997 to allow for a special version to be used to check
	whether there will be some columns with only gaps -- this is used by 
	setting the file parameter to NULL -- no printing is done.
*/
{
   unsigned int* stat;
   int length;
   int i,j,k,l;
	int start;
	unsigned int filter;
	t_sequence *seq;
	int i1,i2;
	int ii;
	t_block* block;
	char match_string[1000];
	int len;


	if (nr==1) {
		Special_Print_Pattern_Matches_For_One_Sequence(file,block_pointer,seqs,nr,Pat_Info);
		return;
	}


	block= (t_block*)block_pointer;
 
   stat= Pat_Info->bit_vector;
   length= Pat_Info->length;
 
	if (file) 
		fprintf (file,"Occurrences: %d(%d)\n",Pat_Info->nr_occ,Pat_Info->nr_match);
 
   for (i=0;i<nr;i++) {
		start= block->start[i];
      for (j=0;j<block->length[i];j++) {
         for (k=0;k<PER_WORD;k++) {
				filter= (1<<k);
            if (stat[start+j] & filter ) {
					for (l=0;l<Pat_Info->num_sub;l++) {
						if (Pat_Info->bit_vecs[l][start+j] & filter) {

							seq= seqs[i];
							i1= (int)(j*PER_WORD+k); 
							len= Match_To_String (Pat_Info,i1,l,match_string,seq,(file!=NULL));
							i2= i1+len-1;
							if (file) {
   							fprintf (file,"%14s%c:%6d-%6d: ",seq->filename,(seq->restricted?'*':' '),i1+1,i2+1);
								for (ii= i1-5;ii<i1;ii++)
      					 		fprintf (file,"%c",(ii>=0 ? (tolower(seq_element(seq,ii))) : ' ') );

								fprintf (file," %s ",match_string);
   							for (ii=i2+1;ii<=i2+5;ii++)
      							fprintf (file,"%c",(ii<seq->length ? (tolower(seq_element(seq,ii))) : ' ') );
   							fprintf(file,"\n");
							}
						}
					}
				}
         }
      }
	}

}


/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
