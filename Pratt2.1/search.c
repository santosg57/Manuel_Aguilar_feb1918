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
#include <time.h>

#include "menu.h"
#include "hit.h"
#include "sequence.h"
#include "al.h"
#include "pattern.h"
#include "block.h"
#include "scan.h"
#define SEARCH_MAIN
#include "search.h"
#include "swiss.h"

#include "tree.h"
#include "mst.h"
#include "divtest.h"

int Graph_Nr_Nodes;

extern char version[];
extern time_t start_sec;
extern int sum_length_sequences;
extern t_sequence* query_sequence;
extern float apr_single[26];
extern float Hit_Needed;

time_t start_search_time;

int Pat_Info_Remove_Gappy_Columns (t_Pat_Info *, t_block* );

/*
#define B_Options->Swiss_Flat_File "/home/tmp/inge/swiss31.dat"
#define B_Options->Swiss_Flat_File "/u2/embnet/mpsrch/dbases/pro/swiss31.dat"
*/

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

/*
#define DEBUG
*/

#ifdef DEBUG
#define CHECK(a) assert(a)
#else
#define CHECK(a)
#endif


#define GAP_RED 0.5

/*
#define GAP_LOG
*/


/**************** VARIABLES GLOBAL TO THIS FILE *************/
static t_Pattern *Patterns_Table[G_MAX_COMP+10];
static unsigned int ****Bit_Vectors;
static unsigned int ***Seq_Hits;
static int Length_Seq_Hits;
static float Threshold_Info;
int Nr_Patterns_Evaluated= 0;
static int Incl_Tab[1000];
static t_Pat_Info* Pat_Being_Refined=NULL;
static unsigned int Non_Flex_Mask;
static int B_Quotient;
static int Nr_Pat_Subtree;
static t_Pat_Node *Start_Node;
static int Pull_Out;

static int Truncated= 0;

/**************** GLOBAL VARIABLES DEFINED HERE *************/
int Length_Block;
int *Nr_Hits;
t_Block_Options *B_Options;
int Nr_Entries_Swiss_Prot;
t_tree *Tree;
float Dist_Matrix[MAX_NUM_DIST][MAX_NUM_DIST];

/**************** GLOBAL VARIABLES DEFINED ELSEWHERE ********/
extern t_sequence** seqs;
extern int nrseqs;
extern int virt_nrseqs;
extern int Min_Nr_Occs;
extern unsigned int Member_In_Set[];
extern int Number_In_Set[];
extern float Freq_Set[];
extern float Bit_Set[];
extern int Nr_Sets;
extern int Nr_Sets_Block;
extern int Nr_Sets_Unit;
extern int Nr_Sets_In_Search;
extern float Shannon;
extern char AA_Tab[];
extern int Length_AA_Set;
extern unsigned int AA_In_Set[26][10];
extern int Num_Ones[];
extern t_Generalize_Symbol_List **Gen_Symbol_Lists;


void flush_header();
void flush_status(t_Hits* );

/**************** LOCAL ROUTINES ****************************/
static int Refine_Pattern(t_block* ,t_Pat_Info*,t_Hits* );
static int Rec_Greedy_Refine_Pattern(t_block* , t_Pattern* ,t_Hits *,t_Search *,int,int );

static int Rec_Find_Pattern_Using_Graph(t_block* , t_Pattern* ,t_Hits *,t_Pat_Graph *,t_Pat_Edge* ,float*,int );
static int Calculate_Cardinality(t_Pattern *,int);
static void Hits_Insert_Pat_Info_Special(t_block* ,t_Hits *,t_Pat_Info *);


/**************** LOCAL MACROES  ****************************/
#define LOOK_MASK 65535
#define NUM_ONES32(n) (Num_Ones[(n)&LOOK_MASK]+Num_Ones[(n)>>16])

t_Hits* Refine_Hits(t_block *block,t_Hits* hits)
/******************************************************************/
{
	t_Hits *Hits2;
	t_Hit_Entry *Entry;
	int i=1;

	printf ("\n\n");

	if (hits->First->Pat==NULL)
		return hits;

	printf ("Refinement phase:\n");

	if (B_Options->Diagnostic)
		Hits2= Init_Hit_List(B_Options->Length_Hit_List,B_Options->Min_Information);
	else
		Hits2= Init_Hit_List(B_Options->Length_Hit_List,B_Options->Min_Information);

	for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next)
	{
		if (!Entry->Pat)
			continue;
		printf ("\rRefining pattern %d   ",i);
		fflush(stdout);
		Refine_Pattern(block,Entry->Pat,Hits2);

		i++;
	}
	return Hits2;
}


static double Gap_Penalty(int i)
{
#ifdef GAP_LOG
	return GAP_RED*(log((double)i+1.0))/log(2.0);
#else
	return GAP_RED*(double)i;
#endif

}

static int Allocate_Structures()
/******************************************************************/
{
	int h,i,j;
	t_Pattern* Pat;
	static int Already_Allocated= 0;
	int num_comp;
	size_t size_comps;

	printf ("  Memory is allocated for the search\n");
	Non_Flex_Mask= 1<<16;

	num_comp= B_Options->Max_Num_Comp;

	size_comps= num_comp*sizeof(int);

	if (Already_Allocated) {
		return 1;
	}

	Already_Allocated= 1;
	
	Nr_Hits= calloc(20,sizeof(int));
	assert (Nr_Hits);

	for (i=0;i<B_Options->Max_Num_Comp;i++)
	{
		Patterns_Table[i]= Pat= malloc(sizeof(t_Pattern));
		assert (Pat);
		Pat->debug= PDEBUG;
		Pat->num_comp= i+1;
		Pat->positions= malloc(size_comps);
		Pat->spec_pos= malloc(size_comps);
		Pat->gaps= malloc(size_comps);
		Pat->flexi= malloc(size_comps);
		assert (Pat->positions && Pat->spec_pos && Pat->gaps && Pat->flexi);
		Pat->num_sub= 0;

		Pat->lengths= malloc(B_Options->Max_Flex_Prod*sizeof(int));
		Pat->bit_vecs= malloc(B_Options->Max_Flex_Prod*sizeof(unsigned int*));
		assert (Pat->lengths && Pat->bit_vecs);


		for (j=0;j<B_Options->Max_Flex_Prod;j++) {
			Pat->bit_vecs[j]= malloc(Length_Block*sizeof(unsigned int));
			assert (Pat->bit_vecs[j]);
		}
	}
	

	Bit_Vectors= calloc(B_Options->Max_Num_Comp,sizeof(unsigned int***));
	assert (Bit_Vectors);


	for (h=0;h<B_Options->Max_Num_Comp;h++)	
	{
		Bit_Vectors[h]= calloc(B_Options->Max_Flex_Prod,sizeof(unsigned int**));
		assert (Bit_Vectors[h]);
		for (i=0;i<B_Options->Max_Flex_Prod;i++)
		{
			unsigned int** vec;
			vec= Bit_Vectors[h][i]= malloc(max(Nr_Sets_Block,B_Options->Max_Gap+1)*sizeof(unsigned int*));
			assert (vec);
			for (j=0;j<=max(Nr_Sets_Block-1,B_Options->Max_Gap);j++) {
				vec[j]= malloc(Length_Block*sizeof(unsigned int));
				assert (vec[j]);
			}
		}
	}

	Length_Seq_Hits= 1+(virt_nrseqs/PER_WORD);
	Seq_Hits= malloc(B_Options->Max_Num_Comp*sizeof(unsigned int**));
	assert (Seq_Hits);

	for (h=0;h<B_Options->Max_Num_Comp;h++)	
	{
		unsigned int** s;
		s= Seq_Hits[h]= 
			malloc(max(Nr_Sets_Block+1,B_Options->Max_Gap+1)*sizeof(unsigned int*));
		assert (s);

		for (i=0;i<=max(B_Options->Max_Gap,Nr_Sets_Block);i++) {
			s[i]= calloc(Length_Seq_Hits,sizeof(unsigned int));
			assert (s[i]);
		}
	}

	return 1;
}


static int Refine_Pattern(t_block* block,t_Pat_Info *Pat_Info,t_Hits* Hits2)
/******************************************************************/
{
	unsigned int* bit_vector;
	t_Pattern *Pat;
	int j,k;
	int success;
	t_Search Search;
	int pos1;
	int nr;

	Pull_Out= 0;   /* set this to one once a maximum pattern is found */

	Pat_Being_Refined= Pat_Info;

	/* Added to avoid trying to refine a pattern that cannot possibly be
		improved   Febr 7th 1996 Inge Jonassen                            */
	if (Pat_Info->num_comp==B_Options->Max_Num_Comp) {
		Hits_Insert_Pat_Info_Special(block,Hits2,Pat_Info);
		return 1;
	}

	Search.Pat= Pat_Info;
	Search.pos= 0;
	Search.gaps= 0;
	Search.flexi= 0;
	
	success= 0;
	Threshold_Info= 1.0;

	pos1= Pat_Info->positions[0];

	Pat= Patterns_Table[0];
	CHECK (Pat->debug==PDEBUG);
	Pat->length= 1;
	Pat->lengths[0]= 1;
	Pat->num_comp= 1;

	CHECK (Pat->positions && Pat->gaps && Pat->flexi);

	Pat->gaps[0]= Pat->flexi[0]= 0; 
	Pat->num_sub= 1;
	Pat->num_flex= 0;
	Pat->num_sub_new= Pat_Info->num_sub;

/*
   for (j=0;j<26;j++)
      if (Member_In_Set[pos1]&(1<<j))
          fprintf (stderr,"%c",(char)(j+'A'));
   fprintf (stderr,"\n");
*/

	Pat->info= Bit_Set[pos1];

	Pat->info_MDL= Pat->info;

	Pat->positions[0]= pos1;
	Pat->spec_pos[0]= Member_In_Set[pos1];

	for (j=0;j<Pat->num_sub_new;j++) {
		for (k=0;k<Length_Block;k++) {
			Pat->bit_vecs[j][k]= Pat_Info->bit_vecs[j][k];
		}
	}
	Pat->cardinality= Calculate_Cardinality(Pat,1);

	Nr_Sets_In_Search= Nr_Sets_Block;

	nr= Hits2->Num;
	success= Rec_Greedy_Refine_Pattern(block,Pat,Hits2,&Search,1,B_Options->Min_Nr_Seqs_Matching);

	if (nr==Hits2->Num)
		Hits_Insert_Pat_Info_Special(block,Hits2,Pat_Info);

	return success;
}



/*
t_Hits *Search_For_Patterns(t_block* block)
{
	t_Pattern *Pat;
	int i,j;
	int success;
	t_Hits* Hits;

	success= 0;

	Hits= Init_Hit_List(B_Options->Length_Hit_List,B_Options->Min_Information);

	Pat= Patterns_Table[0];
	CHECK (Pat->debug==PDEBUG);
	Pat->length= 1;
	Pat->lengths[0]= 1;
	Pat->num_comp= 1;
	assert (Pat->positions && Pat->gaps && Pat->flexi);
	Pat->gaps[0]= Pat->flexi[0]= 0; 
	Pat->num_sub= 1;
	Pat->num_flex= 0;

	for (i=0;i<Nr_Sets_In_Search;i++)
	{
		CHECK (Bit_Set[i]>1.0);
	
		if (Bit_Set[i]<Threshold_Info)
			continue;

      for (j=0;j<26;j++)
         if (Member_In_Set[i]&(1<<j))
            fprintf (stderr,"%c",(char)(j+'A'));
      fprintf (stderr,"\n");

		Pat->info= Bit_Set[i];
		Pat->info_MDL= Pat->info;
		Pat->positions[0]= i;
		Pat->spec_pos[0]= Member_In_Set[i];
		Pat->gaps[0]= 0;
		for (j=0;j<Length_Block;j++) {
			Pat->bit_vecs[0][j]=  block->vectors[i][j];
		}
		success= Rec_Guaranteed_Find_Pattern(block,Pat,Hits,NULL,0) || success;
	}
	return Hits;
}
*/


int nr_total,nr_now;


t_Hits *Find_Conserved_Patterns_Using_Graph(t_block* block, t_Pat_Graph* Graph)
/******************************************************************/
{
	t_Pattern *Pat;
	int i,j,k,l;
	int success;
	t_Hits* Hits;
	t_Pat_Node *Node;
	t_Pat_Edge *Edge,*Edge_Out;
	int symbol;
	float max_info;
	float max_score;
	unsigned int symbol_set;
	int symbol_index;
	char buf[20];
	unsigned int *start_vector;
	int used_time;
	time_t now;

	success= 0;
	Graph_Nr_Nodes= Graph->nr_nodes;

	Hits= Init_Hit_List(B_Options->Length_Hit_List,B_Options->Min_Information);

	Pat= Patterns_Table[0];
	CHECK (Pat->debug==PDEBUG);
	Pat->length= 1;
	Pat->lengths[0]= 1;
	Pat->num_comp= 1;
	assert (Pat->positions && Pat->gaps && Pat->flexi);
	Pat->gaps[0]= Pat->flexi[0]= 0; 
	Pat->num_sub= 1;
	Pat->num_flex= 0;
	Pat->num_seq= virt_nrseqs;

	Edge= malloc(sizeof(t_Pat_Edge));

	if (B_Options->Restrictions_Input==start) {
		start_vector= calloc(block->len,sizeof(unsigned int));
		for (i=0;i<nrseqs;i++)
			if (seqs[i]->restricted) {
				for (j=0;j<seqs[i]->restricted;j++)
					for (k=seqs[i]->first[j];k<=seqs[j]->last[0];k++)
						start_vector[block->start[i]+(k/PER_WORD)]|= 1<<(k%PER_WORD);
			}
			else {
				for (k=0;k<seqs[i]->length;k++)
					start_vector[block->start[i]+(k/PER_WORD)]|= 1<<(k%PER_WORD);
			}
	}
		

/*
	for (i=0;i<Graph->nr_nodes;i++)
*/

	flush_header();


	nr_total= Graph->nr_nodes; nr_now= 0;
	for (i=Graph->nr_nodes-1;i>=0;i--)
	{
		max_info= 0.0;
		nr_now++;
/*
		if (Graph->nodes[i]->seq_index>(nrseqs-B_Options->Min_Nr_Seqs_Matching))
			continue;
*/

		B_Quotient= B_Options->Quotient;
		Nr_Pat_Subtree= 0;
		Node= Graph->nodes[i];
		Start_Node= Node;
		success= 0;

		Node->Analysed= 0;
		Node->Best_Pat= NULL;

		now= time(NULL);
		used_time= (now-start_sec);
		if (used_time>=B_Options->Max_Time) {
			Truncated= 1;
			return Hits;
		}


/*
		if (Node->length_longest_path<3)
			continue;
*/

		if ( Node->Max_Score[B_Options->Max_Num_Comp-1][B_Options->Max_Num_Flex]
			 < B_Options->Min_Information )
		{
			continue;
		}

		symbol_set= Member_In_Set[Node->symbol];

/*
     	printf ("\n%4d/%4d %8.2f\n",i,Graph->nr_nodes,Node->Max_Score[B_Options->Max_Num_Comp-1][B_Options->Max_Num_Flex]);
*/

		for (symbol_index=0;
			  symbol_index<Gen_Symbol_Lists[Node->symbol]->nr_symbols;
			  symbol_index++)
		{
			symbol= Gen_Symbol_Lists[Node->symbol]->symbols[symbol_index];

			assert ((symbol_set&Member_In_Set[symbol])==symbol_set);

			/* Initialize pattern starting in Node */

			Pat->info= Bit_Set[symbol];

			if (B_Options->MDL_flag)
				Pat->info_MDL= Pat->info-B_Options->MDL_constant0- ((B_Options->MDL_constant1*Pattern_Num_Char(Pat)+B_Options->MDL_constant2*Pattern_Num_X(Pat)+B_Options->MDL_constant3)/((float)max(1,Pat->num_seq)));
			else
				Pat->info_MDL= Pat->info;

			Pat->positions[0]= symbol;
			Pat->spec_pos[0]= Member_In_Set[symbol];
			Pat->gaps[0]= 0;
			Pat->num_sub= 1;
			Pat->num_flex= 0;

			if (B_Options->Restrictions_Input==start) {
				for (j=0;j<block->len;j++)
					Pat->bit_vecs[0][j]= start_vector[j] & block->vectors[symbol][j];
			}
			else
				memcpy(Pat->bit_vecs[0],block->vectors[symbol],Length_Block*sizeof(unsigned int));

			Pat->cardinality= Calculate_Cardinality(Pat,0);

/*
			for ( max_info= 0.0, j=0;j<Node->nr_edges;j++) 
*/
			for ( j=0;j<Node->nr_edges;j++) 
			{
				int flex;
				int flex_out;

				Edge_Out= Node->edges[j];

				flex= Edge_Out->max_aa-Edge_Out->min_aa;
				flex_out= (flex>0 ? 1 : 0);

				max_score= Pat->info + Edge_Out->Node->Max_Score[B_Options->Max_Num_Comp-2][B_Options->Max_Num_Flex-flex_out] - Gap_Penalty(flex);

/*
				if (Edge_Out->Node->Analysed) {
					max_score=min(max_score,Pat->info + Edge_Out->Node->Max_Score_Obtained - Gap_Penalty(flex));
				}
*/
	
				if (max_score <= max_info || max_score < B_Options->Min_Information)
				{
					continue;
				}
	
				success= Rec_Find_Pattern_Using_Graph(block,Pat,Hits,Graph,Edge_Out,&max_info,symbol_index) || success;
				if (success) {
					Node->Analysed= 1;
					Node->Max_Score_Obtained= max_info;
				}
			}
		}
	}

	if (B_Options->Restrictions_Input==start) 
		free (start_vector);

	flush_status(Hits);

	printf ("\n");
	
	return Hits;
}

static void Score_Pat_Info(t_Pat_Info *pat,t_block* block)
/******************************************************************/
/* This is to be used by Remove_Gappy_Columns to set the score right
	after the unnecassary flexibility has been removed
	Added 4/2/97 IJ
*/
{
	int i;

	for (i=0;i<pat->num_comp;i++) {
		if (pat->nflexi[i]<pat->flexi[i])
			pat->information= pat->information + Gap_Penalty(pat->flexi[i]) - Gap_Penalty(pat->nflexi[i]);
	}

/*
	if (B_Options->MDL_flag)  {
		pat->fitness= pat->info_MDL= pat->information- B_Options->MDL_constant0-((B_Options->MDL_constant1*Pat_Info_Num_Char(pat)+B_Options->MDL_constant2*Pat_Info_Num_X(pat)+B_Options->MDL_constant3)/((float)max(1,pat->num_seq)));
	}
	else if (!Pat_Info_Divergence(pat,block))
		pat->fitness= pat->information;
*/
}

static void Hits_Insert(t_block* block,t_Hits *hits,t_Pattern *Pattern,int initialsearch_flag)
/******************************************************************/
{
	t_Pat_Info *pat;

	/* if the conditions are fullfilled this pattern will not be inserted
		and there is no need to make a Pat_Info structure from it:        */
	if ((B_Options->MDL_flag==0) && 
		 (B_Options->Tree_Input==0) && 
		 (B_Options->Dist_Input==0) && 
		 !Hit_List_Will_Insert_Entry(hits,Pattern->info))
		return;

	if (initialsearch_flag) {
		if (!Pattern_Already_In_List(hits,Pattern)) {
			pat= Pat_Info_From_Pattern(Pattern,initialsearch_flag);

			Pat_Info_Remove_Gappy_Columns(pat,block);
			
			if (B_Options->MDL_flag) {
				pat->info_MDL= pat->information- B_Options->MDL_constant0-((B_Options->MDL_constant1*Pat_Info_Num_Char(pat)+B_Options->MDL_constant2*Pat_Info_Num_X(pat)+B_Options->MDL_constant3)/((float)max(1,Pattern->num_seq)));
				pat->fitness= pat->info_MDL;
			}
			else if (!Pat_Info_Divergence(pat,block))
				pat->fitness= pat->information;

			if (Hit_List_Will_Insert_Entry(hits,pat->fitness)) {
				Hit_List_Insert_Entry(hits,pat,pat->fitness,0);
			}
			/*
				this caused the program to crash at a linux machine
			else
				Pat_Info_Free(pat);
			*/
		}
	}
	else
	{
		if (Pattern->num_comp==B_Options->Max_Num_Comp)
			Pull_Out= 1;

		pat= Pat_Info_From_Pattern(Pattern,initialsearch_flag);

		Pat_Info_Remove_Gappy_Columns(pat,block);

		pat->refined_from= Pat_Being_Refined;
		if (B_Options->Diagnostic) {
			Scan_Swiss_Prot_Refined_Pattern(pat);
		}
		else {
			pat->nr_occ_swiss= 0;
			pat->nr_match_swiss= 0;
		}
		Statistics_Refined_Pattern(block,pat);
		pat->num_seq= pat->nr_match;
		Hit_List_Insert_Entry(hits,pat,pat->fitness,1);
	}
}


static void Hits_Insert_Pat_Info_Special(t_block* block,t_Hits *hits,t_Pat_Info *pat)
/******************************************************************/
/* Used for inserting directly a pattern (to be refined) already  */
/* having the max. number of components allowed.                  */
/* To avoid waisting time trying to refine a "perfect thing"      */
/* Added Febr 7th 1996  Inge Jonassen                             */
/******************************************************************/
{
	pat->refined_from= Pat_Being_Refined;
	if (B_Options->Diagnostic) {
		Scan_Swiss_Prot_Refined_Pattern(pat);
	}
	else {
		pat->nr_occ_swiss= 0;
		pat->nr_match_swiss= 0;
	}
	Statistics_Refined_Pattern(block,pat);
	Hit_List_Insert_Entry(hits,pat,pat->fitness,1);
}



static int More_Occ(const void *i1, const void *i2)
/******************************************************************/
{
	return (Nr_Hits[(*(int*)i2)]-Nr_Hits[(*(int*)i1)]);
}


static int Nr_Seqs_Hit(t_block *block, unsigned int *vector)
/******************************************************************/
{
	int j, num,prev;

	prev= (-1);
  	num= 0;
  	for (j=0;j<Length_Block;j++) {
		if (vector[j] && (block->seq[j]!=prev)) {
     		num++;
     		prev= block->seq[j];
  		}
	}

	return num;
}



static int Rec_Greedy_Refine_Pattern(t_block* block, t_Pattern* Pat,t_Hits *hits,t_Search *Search,int move,int min_seq_match)
/******************************************************************/
{
	unsigned int* vec1;
	unsigned int* end1;
	unsigned int* vec2;
	unsigned int* end2;
	unsigned int* vec3;
	int *ipt1, *ipt2;
	int j,k,l,n;
	int gap; 
	t_Pattern *Local;
	int success;
	int symbol;
	int sub;
	int seq_index;
	int max_gap;
	unsigned int ***Bit_Vec;
	unsigned int **Seq_Hit;
	int search_pos, search_gaps, search_flexi;
	static int *order_aa=NULL;
	int nr_hit;
	float info;
	unsigned int union_seqs;
	int num, prev;
	int nr_included;
	int ok;
	static unsigned int Set_Vector[10];
	int *intvec1, *intend1;
	size_t length_bit_vector;
	int tt;

	if (Pull_Out)
		return 0;

	if (order_aa==NULL)
		order_aa= (int*)calloc(20,sizeof(int));

	assert(order_aa);

	assert (Search);

	search_pos= Search->pos;
	search_gaps= Search->gaps;
	search_flexi= Search->flexi;

	if (move) {
		Search->pos++;
		Search->gaps= 0;
		Search->flexi= 0;
	}


	success= 0;


	if (Pat->num_comp>=(B_Options->Max_Num_Comp) || (Pat->length>=(B_Options->Max_Length)))
		return 0;

	Local= Patterns_Table[Pat->num_comp];

	Bit_Vec= Bit_Vectors[Pat->num_comp-1];
	Seq_Hit= Seq_Hits[Pat->num_comp-1];

	if (!Pat_Copy(Pat,Local)) {
		printf ("Error: Cannot copy Pattern %d\n",Pat->num_comp);
		exit (1);
	}
	length_bit_vector= Length_Block*sizeof(unsigned int);

	if ((Search->pos < Search->Pat->num_comp) && (Search->Pat->flexi[Search->pos]>0) )
	{
		/* Flexible wildcard -- cannot refine this */
		assert (Search->flexi==0);
		assert (Search->gaps==0);
		tt= Local->positions[Local->num_comp]= Search->Pat->positions[Search->pos];
		Local->spec_pos[Local->num_comp]= Member_In_Set[tt];
		Local->gaps[Local->num_comp]= Search->Pat->gaps[Search->pos]-Search->gaps;
		Local->flexi[Local->num_comp]= Search->Pat->flexi[Search->pos]-Search->flexi;
		Local->num_sub*= (1+Search->Pat->flexi[Search->pos]);
		Local->length= Pat->length + Local->gaps[Local->num_comp] +
                                   Local->flexi[Local->num_comp] + 1;
		for (j=0;j<Local->num_sub_new;j++) {
			if (j<Local->num_sub)
				Local->lengths[j]= Pat->lengths[j%Pat->num_sub] + 
                               Local->gaps[Local->num_comp]+ j/Pat->num_sub + 1;
			else
				Local->lengths[j]= Local->lengths[j%Local->num_sub];

         vec1= Pat->bit_vecs[j];
         vec2= Local->bit_vecs[j];
			memcpy(vec2,vec1,length_bit_vector);
		}
		
		Local->info= Pat->info + 
                 	Bit_Set[Local->positions[Local->num_comp]] - 
						Gap_Penalty(Local->flexi[Local->num_comp]);

		Local->num_comp++;

		if (!Rec_Greedy_Refine_Pattern(block,Local,hits,Search,1,min_seq_match))
			Hits_Insert(block,hits,Local,0);

		return 1;
	}

	if (Search->pos >= Search->Pat->num_comp)
		max_gap= min(B_Options->Max_Gap,B_Options->Max_Length-(Pat->length+1));
	else	
		max_gap= Search->Pat->gaps[Search->pos]-Search->gaps;

	for (gap=0;gap<max_gap;gap++)
	{
		for (j=0;j<Nr_Sets_Unit;j++)
			memset(Seq_Hit[j],0,Length_Seq_Hits*sizeof(unsigned int));

		for (symbol=0;symbol<Nr_Sets_Unit;symbol++)
		{
			for (sub=0;sub<Pat->num_sub_new;sub++) {
           	vec1= block->vectors[Nr_Sets_Block*(Pat->lengths[sub%(Pat->num_sub)]+gap)+symbol];
				vec2= Pat->bit_vecs[sub];
				vec3= Bit_Vec[sub][symbol];
				end1= vec1+Length_Block;
				while (vec1<end1)
					*vec3++= *vec1++ & *vec2++;
			}
			k= 0;
			while (k<Length_Block)
			{
				sub= 0;
				while ((sub<Pat->num_sub_new) && (Bit_Vec[sub][symbol][k]==0))
					sub++; 
				if (sub<Pat->num_sub_new)
				{
					CHECK (Bit_Vec[sub][symbol][k]!=0);
					seq_index= block->seq[k];
					CHECK(seq_index>=0 && seq_index<virt_nrseqs);
					Seq_Hit[symbol][seq_index/PER_WORD]|= (1<<(seq_index%PER_WORD));
					k= block->start[seq_index+1];
				}
				else k++;
			}
		}
		assert (Nr_Sets_Unit<=20);

		for (j=0;j<Nr_Sets_Unit;j++)
		{
			for (Nr_Hits[j]=k=0;k<Length_Seq_Hits;k++) {
				Nr_Hits[j]+= Num_Ones[Seq_Hit[j][k] & LOOK_MASK]+Num_Ones[Seq_Hit[j][k]>>16];
			}
			CHECK(Nr_Hits[j]<=virt_nrseqs);
		}

		for (j=0;j<Nr_Sets_Unit;j++)
			order_aa[j]= j;

		CHECK(Nr_Sets_Unit<=20);

		qsort((void*)order_aa,(size_t)Nr_Sets_Unit,sizeof(int),More_Occ);

#ifdef DEBUG
	for (k=0;k<20;k++)
		assert (order_aa[k]<=20);
#endif

		memset(Seq_Hit[Nr_Sets_Unit],0,Length_Seq_Hits*sizeof(unsigned int));

		nr_included= nr_hit= 0;
		memcpy(Set_Vector,AA_In_Set[AA_Tab[order_aa[0]]-'A'],Length_AA_Set*sizeof(unsigned int));

		j= 0;
		while (nr_hit<min_seq_match && j<Nr_Sets_Unit)
		{
			ok= 0;
			vec1= Set_Vector;
			vec2= AA_In_Set[AA_Tab[order_aa[j]]-'A'];
			end1= vec1+Length_AA_Set;
			while (vec1<end1)
				if (*vec1++ & *vec2++)
					ok= 1;

			if (ok)
			{
				Incl_Tab[nr_included]= order_aa[j];
				nr_included++;

				vec1= Set_Vector;
				vec2= AA_In_Set[AA_Tab[order_aa[j]]-'A'];
				while (vec1<end1)
					*vec1++&= *vec2++;

				nr_hit= 0;
				vec1= Seq_Hit[Nr_Sets_Unit];
				vec2= Seq_Hit[order_aa[j]];
				end1= vec1+Length_Seq_Hits;
				while (vec1<end1) {
					*vec1 |= *vec2;
					nr_hit+= Num_Ones[(*vec1) & LOOK_MASK]+ Num_Ones[(*vec1) >> 16];
					vec1++; vec2++;
				}
			}
			j++;	
			CHECK(nr_hit <= virt_nrseqs);
		}

		if (nr_hit>=min_seq_match && nr_included>=1)
		{
			float freq;
			int ii,jj,kk;
			int keepon;

			/* First find the amino acid set with the most info content that contains all
				the amino-acids seen in this column */
			float best_info;
			int best_index;

			best_info= 0.0; best_index= (-1);
			for (ii=0;ii<Length_AA_Set;ii++) 
			{
				if (Set_Vector[ii]) {
					for (jj=0;jj<PER_WORD;jj++) {
						kk=ii*PER_WORD+jj;
						if (Set_Vector[ii]&(1<<jj)) {
							if (Bit_Set[kk]>best_info) {
								best_info= Bit_Set[kk];
								best_index= kk;
							}
						}
					}
				}
			}
			assert(best_index>=0);

			if (B_Options->Refine_Generalise) {
				/* Include the ones in this set that has not already been included */
				for (ii=0;ii<Nr_Sets_Unit;ii++) {
					if (Member_In_Set[best_index]&(1<<AA_Tab[ii]-'A')) {
						for (ok=jj=0;jj<nr_included;jj++) {
							if (Incl_Tab[jj]==ii)
								ok=1;
						}
						if (!ok)
							Incl_Tab[nr_included++]=ii;
					}
				}
			}

			/* calculate info -- the info content of the new added pattern component:
				this should be the same as best_info when a complete set from the allowed ones
				(in Pratt.sets) is included */

			freq= 0.0;
			intvec1= Incl_Tab;
			intend1= Incl_Tab+nr_included;
			while (intvec1< intend1)
				freq+= Freq_Set[*intvec1++];

			info= Shannon;

			intvec1= Incl_Tab;
			intend1= Incl_Tab+nr_included;
			while (intvec1<intend1)
			{
				float freq_j;
				assert (*intvec1<Nr_Sets_Unit);
				freq_j= Freq_Set[*intvec1++]/freq;
				info+= freq_j*log(freq_j)/log(2.0);
			}

			n= Pat->num_comp;
			assert (n>0);
			Local->num_comp= n+1;
			Local->positions[n]= -1;
			vec1= &(Local->spec_pos[n]);
			*vec1= 0;
			ipt1= Incl_Tab; ipt2= Incl_Tab+nr_included;
			while (ipt1<ipt2)
				*vec1 |=  Member_In_Set[*ipt1++];
				/*
				*vec1 |= 1 << *ipt1++;
				*/

			Local->gaps[n]= gap;
			Local->flexi[n]= 0;

			Local->length= Pat->length + gap + 1;
			Search->gaps+= gap+1;

			Local->info= Pat->info+info;
			Local->num_sub= Pat->num_sub;

			for (j=0;j<Local->num_sub_new;j++)
			{
				unsigned int **temp_pt;
				CHECK (n+1==Local->num_comp);
				Local->lengths[j]= Pat->lengths[j%Pat->num_sub]+gap+1;
				vec1= Local->bit_vecs[j];
				memset(vec1,0,Length_Block*sizeof(unsigned int));
				temp_pt= Bit_Vec[j];
				for (l=0;l<nr_included;l++) {
					vec1= Local->bit_vecs[j];
					end1= vec1+Length_Block;
					vec2= temp_pt[Incl_Tab[l]];
					while (vec1<end1)
						*vec1++ |= *vec2++;
				}
			}

			Local->cardinality= Calculate_Cardinality(Local,1);

			keepon= 0;

			if ((success=Rec_Greedy_Refine_Pattern(block,Local,hits,Search,0,min_seq_match))==0) 
			{ 
				if (Local->num_comp==B_Options->Max_Num_Comp && Search->pos<Search->Pat->num_comp) {
				   /* set to one if max number of components is exhausted and
						complete pattern (to be refined) has not been included.
						In order to force inclusion of complete pattern (to be
						refined).
					*/
					keepon= 1; 
				}
				else 	if (Search->pos>=Search->Pat->num_comp) {
					success= 1;
					Hits_Insert(block,hits,Local,0);
				}
			}

			if (Pull_Out)
				return success;

			/*
			if (!keepon && (B_Options->Quotient!=0) && 
			*/
			if (!keepon && (B_Options->Quotient!=0) && 
				 (Local->cardinality>=(Pat->cardinality/B_Options->Quotient)))
				return success;
			
		}

		Search->gaps= (move?0:search_gaps);
		Search->flexi= (move?0:search_flexi);
	}

	if (Search->pos >= Search->Pat->num_comp)
		return success;

	Local->num_comp= Pat->num_comp;
	tt= Local->positions[Local->num_comp]= Search->Pat->positions[Search->pos];
	Local->spec_pos[Local->num_comp]= Member_In_Set[tt];
	Local->gaps[Local->num_comp]= Search->Pat->gaps[Search->pos] - Search->gaps;
	Local->flexi[Local->num_comp]= 0;
	Local->num_sub= Pat->num_sub;
	Local->num_sub_new= Pat->num_sub_new;
	Local->length= Pat->length + Local->gaps[Local->num_comp]+Local->flexi[Local->num_comp]+1;
	for (j=0;j<Local->num_sub_new;j++) {
		Local->lengths[j]= Pat->lengths[j%Pat->num_sub] + Local->gaps[Local->num_comp] + 1;
		memcpy(Local->bit_vecs[j],Pat->bit_vecs[j],Length_Block*sizeof(unsigned int));
	}

	Local->info= Pat->info + Bit_Set[Local->positions[Local->num_comp]];
	Local->num_comp++;

	if ((success=Rec_Greedy_Refine_Pattern(block,Local,hits,Search,1,min_seq_match)==0) &&
		 (Search->pos>=Search->Pat->num_comp)) {
		success= 1;
		Hits_Insert(block,hits,Local,0);
	}

#ifdef DEBUG
	Pat_Check(Pat);
	Pat_Check(Local);
#endif

	Search->pos= search_pos;
	Search->gaps= search_gaps;
	Search->flexi= search_flexi;

	return success;
}




static int Rec_Find_Pattern_Using_Graph(t_block* block, t_Pattern* Pat,t_Hits *hits,t_Pat_Graph *Graph,t_Pat_Edge* Edge,float *Max_Score_Seen,int last_symbol)
/******************************************************************/
{
	unsigned int* vec1;
	unsigned int* vec2;
	unsigned int* vec3;
	unsigned int* end1;
	int i,j,jj,k,l,m,n;
	int gap; 
	t_Pattern *Local;
	int num_seq;
	int success;
	int symbol;
	int sub;
	unsigned int seq_mask;
	int seq_index;
	int flex;
	int min_dist, max_dist,max_flex;
	int pat_j;
	unsigned int ***Bit_Vec;
	unsigned int **Seq_Hit;
	int num_sub;
	t_Pat_Node* Node;
	t_Pat_Edge* Edge_Out;
	unsigned int symbol_set;
	int found;
	int current_gap;
	float max_score;
	int symbol_index;
	int max_cardinality;
/*
	t_Pat_Lead_Info *Best_Pat;
*/

	found= 0;
	max_cardinality= 0;

	Node= Edge->Node;
	gap= Edge->min_aa;
	flex= Edge->max_aa - Edge->min_aa;


	if (Pat->num_comp==(B_Options->Max_Num_Comp) || (Pat->length>=(B_Options->Max_Length)))
		return 0;
	if (Pat->length+gap+flex>=B_Options->Max_Length)
		return 0;
	if (Pat->num_sub*(flex+1)>B_Options->Max_Flex_Prod)
		return 0;
	if (Pat->num_flex+(flex>0?1:0)>B_Options->Max_Num_Flex)
		return 0;

	Local= Patterns_Table[Pat->num_comp];

	Bit_Vec= Bit_Vectors[Pat->num_comp-1];
	Seq_Hit= Seq_Hits[Pat->num_comp-1];

	if (!Pat_Copy(Pat,Local)) {
		printf ("Error: Cannot copy Pattern %d\n",Pat->num_comp);
		exit (1);
	}

	num_sub= Pat->num_sub;
	symbol_set= Member_In_Set[Node->symbol];

	max_flex= min(B_Options->Max_Flex, (B_Options->Max_Flex_Prod/num_sub)-1);
	max_flex= min(max_flex,B_Options->Max_Gap-Edge->min_aa);
	if (Pat->num_flex>=B_Options->Max_Num_Flex)
		max_flex= 0;

	min_dist= max(0,Edge->max_aa-max_flex);
	max_dist= min(B_Options->Max_Length - Pat->length-1,Edge->min_aa+max_flex);

	for (symbol_index=0;
		  (symbol_index<Gen_Symbol_Lists[Node->symbol]->nr_symbols) &&
		  (B_Options->Quotient==0 || 
			max_cardinality<max(B_Options->Min_Nr_Seqs_Matching,Pat->cardinality/B_Quotient));
		  symbol_index++)
	{
		unsigned int poss;

		if ((poss=Edge->possible[last_symbol][symbol_index])==0)
			continue;

		if (max_flex==0 && (poss&Non_Flex_Mask)==0)
			continue;

		symbol= Gen_Symbol_Lists[Node->symbol]->symbols[symbol_index];

		assert ((symbol_set&Member_In_Set[symbol])==symbol_set);

		for (j=min_dist;j<=max_dist;j++)
			memset(Seq_Hit[j-min_dist],0,Length_Seq_Hits*sizeof(unsigned int));

		for (j=min_dist;j<=max_dist;j++)
		{
			jj= j-min_dist;
			for (sub=0;sub<num_sub;sub++) {
      		vec1= block->vectors[Nr_Sets_Block*(Pat->lengths[sub]+j)+symbol];
				vec2= Pat->bit_vecs[sub];
				vec3= Bit_Vec[sub][jj];
				end1= vec1+Length_Block;
				while (vec1<end1)
					*vec3++= *vec1++ & *vec2++;
			}
			k= 0;
			while (k<Length_Block)
			{
				sub= 0;
				while (sub<num_sub && Bit_Vec[sub][jj][k]==0)
					sub++; 
				if (sub<num_sub)
				{
					seq_index= block->seq[k];
					CHECK(seq_index>=0 && seq_index<virt_nrseqs);
					Seq_Hit[jj][seq_index/PER_WORD]|= (1<<(seq_index%PER_WORD));
					k= block->start[seq_index+1];
				}
				else k++;
			}
		}

		for (k=0;(k<=(max_flex-flex))&&(B_Options->Quotient==0 || max_cardinality<max(B_Options->Min_Nr_Seqs_Matching,Pat->cardinality/B_Quotient));k++) 
		{
			for (l=0;(l<=k) && (B_Options->Quotient==0 || max_cardinality<max(B_Options->Min_Nr_Seqs_Matching,Pat->cardinality/B_Quotient));l++)
			{
				current_gap= gap-l;

				if (current_gap>=0 && current_gap<min_dist)
				{
					printf ("current_gap=%d, min_dist=%d\n",current_gap,min_dist);
					exit (1);
				}

				if (current_gap<0 || current_gap+k+flex>max_dist)
					continue;

				max_score= Pat->info
				  	+ max(0, Edge->Node->Max_Score[B_Options->Max_Num_Comp-Pat->num_comp-1]
					[B_Options->Max_Num_Flex - Pat->num_flex - (flex+k>0 ? 1 : 0)])
				  	- Gap_Penalty(flex+k);

/*
				if (Edge->Node->Analysed) {
					max_score=min(max_score,Pat->info + Edge->Node->Max_Score_Obtained - Gap_Penalty(flex+k));
				}
*/
	

/*
				if ((max_score <= *Max_Score_Seen) || (max_score < B_Options->Min_Information)) {
*/
				if ((max_score <= *Max_Score_Seen) || (max_score < Hit_Needed)) {
					continue;
				}

				for (num_seq=m=0;m<Length_Seq_Hits;m++) {
					seq_mask= 0;
					for (n=0;n<=(flex+k);n++)
						seq_mask|= Seq_Hit[current_gap+n-min_dist][m];
					num_seq+=  NUM_ONES32(seq_mask);
				}
	
				CHECK(num_seq<=virt_nrseqs);
	
				if (num_seq<B_Options->Min_Nr_Seqs_Matching)
				{
					continue;
				}

				found= 1;

				n= Pat->num_comp;
				CHECK (n>0);

				Local->num_seq= num_seq;

				Local->num_comp= n+1;
				Local->positions[n]= symbol;
				Local->spec_pos[n]= Member_In_Set[symbol];
				Local->gaps[n]= current_gap;
				Local->flexi[n]= flex+k;
				Local->length= Pat->length + current_gap + flex+k+ 1;
	
				Local->info= Pat->info+Bit_Set[symbol]-Gap_Penalty(flex+k);

				if (B_Options->MDL_flag)
					Local->info_MDL= Local->info-B_Options->MDL_constant0- ((B_Options->MDL_constant1*(float)Pattern_Num_Char(Pat)+B_Options->MDL_constant2*(float)Pattern_Num_X(Pat)+B_Options->MDL_constant3)/((float)max(1,Local->num_seq)));

				Local->num_sub= Pat->num_sub*(flex+k+1);
				Local->num_flex= Pat->num_flex + ((flex+k)>0 ? 1 : 0);
				for (j=0;j<Local->num_sub;j++)
				{
					pat_j= j%num_sub;
					CHECK (n+1==Local->num_comp);
					Local->lengths[j]= Pat->lengths[pat_j]+current_gap+(j/num_sub)+1;
					vec1= Local->bit_vecs[j];
					vec2= Bit_Vec[pat_j][current_gap-min_dist+(j/num_sub)];
					memcpy(vec1,vec2,Length_Block*sizeof(unsigned int));
				}

				if (Local->info > *Max_Score_Seen) {
					*Max_Score_Seen= Local->info;
/*
					if (Local->info>40.0)
						Start_Node->Best_Pat= Lead_Info_From_Pat(Local,Start_Node->Best_Pat);
*/
				}

				Local->cardinality= Calculate_Cardinality(Local,0);
				max_cardinality= max(max_cardinality,Local->cardinality);
			
				Nr_Patterns_Evaluated++;

				if ((Nr_Patterns_Evaluated%100)==0)
					flush_status(hits);
			
				if (Local->num_comp==B_Options->Max_Num_Comp)
					success= 0;
				else
				{
					for (success=j=0;
					  	j<Node->nr_edges && (B_Options->Quotient==0 || success < max(B_Options->Min_Nr_Seqs_Matching,Local->cardinality/B_Quotient)) ;
					  	j++)
					{
						Edge_Out=Node->edges[j];
	
						max_score= Local->info
					  		+ Edge_Out->Node->Max_Score[B_Options->Max_Num_Comp-Local->num_comp-1]
							[B_Options->Max_Num_Flex - Local->num_flex]
					  		- Gap_Penalty(Edge_Out->max_aa-Edge_Out->min_aa);

/*
						if ((max_score <= *Max_Score_Seen) || (max_score < B_Options->Min_Information))
					;
						else 
*/
							success=Rec_Find_Pattern_Using_Graph(block,Local,hits,Graph,Edge_Out,Max_Score_Seen,symbol_index);

					}
				}
	
				if (!success)
				{
					if (Local->info >= (*Max_Score_Seen-0.01))
						Hits_Insert(block,hits,Local,1);

					if ((Local->info>B_Options->Min_Information) &&
						 (++Nr_Pat_Subtree>=B_Options->Quotient_Inc))
					{
						B_Quotient++;
						Nr_Pat_Subtree= 0;
					}
				}


			} /* l */


		} /* k */

	} /* symbol_index */
	

	return max_cardinality;
}

int Max_Comp_Used(t_Hits *hits)
/******************************************************************/
{
	t_Hit_Entry *temp;
	t_Pat_Info *pat;

	for (temp= hits->First; (temp!=NULL); temp= temp->Next) { 
		if ((pat=temp->Pat) && 
		    (pat->num_comp==B_Options->Max_Num_Comp))
			return 1;
	}

	return 0;
}



int Pat_Info_Remove_Gappy_Columns (t_Pat_Info *pat, t_block* block)
/******************************************************************/
{
	int i,j;
	int fi;
	int flag;

	if (pat->num_sub==1)
		return 0;

	Special_Print_Pattern_Matches(NULL,(void*)block,seqs,nrseqs,pat);

	fi=0;
	flag= 0;
	for (i=0;i<pat->num_comp;i++) {

		if (pat->flexi[i]==0) {
			pat->nflexi[i]= 0;
			continue;
		}
		j= pat->flexi[i];
		while ( (j>=0) && ((pat->used[fi]&(1<<j))==0))
			j--;
		/* 
			Discovered unused flexibility (ie, column with only gaps)
			Remove and updata pattern data.
		*/
		if (j<pat->flexi[i]) {
			pat->special_flag= 1;
			pat->nflexi[i]= j;
			flag= 1;
		}
		fi++;
	}
	if (flag) 
		Score_Pat_Info(pat,block);

	return flag;
}



int Hits_Remove_Gappy_Columns (t_Hits *hits,t_block* block)
/******************************************************************/
/* 
	Inge January 3rd 1997
	This function is to remove not needed flexibility -- this kind
	of flexibility shows in the alignments as columns containing
	only gaps. January 12th -- decided to just give a new flexibility
	array nflexi that contains the needed flexibility. This way I dont
	have to move the bitvectors around. See Special_Print_Pattern_Matches
	in pattern.c for how this is handled.
*/
{
	t_Hit_Entry *temp;
	t_Pat_Info *pat;
	int i,j;
	int fi;
	int flag;

	for (temp= hits->First; (temp!=NULL); temp= temp->Next) 
	{ 
		if ((pat=temp->Pat)==NULL)
			continue;
		if (pat->num_sub==1)
			continue;

		Special_Print_Pattern_Matches(NULL,(void*)block,seqs,nrseqs,pat);

		fi=0;
		flag= 0;
		for (i=0;i<pat->num_comp;i++) {

			if (pat->flexi[i]==0) {
				pat->nflexi[i]= 0;
				continue;
			}
			j= pat->flexi[i];
			while ( (j>=0) && ((pat->used[fi]&(1<<j))==0))
				j--;
			/* 
				Discovered unused flexibility (ie, column with only gaps)
				Remove and updata pattern data.
			*/
			if (j<pat->flexi[i]) {
				pat->special_flag= 1;
				pat->nflexi[i]= j;
				flag= 1;
			}
			fi++;
		}
		if (flag) 
			Score_Pat_Info(pat,block);
	}
}



void Find_Motif_Block(t_Block_Options *Options)
/******************************************************************/
{	
	t_block *block;
	t_Hits *Hits1;
	t_Hits *Hits2;
	t_alignment* al;
	t_Pat_Graph* Pat_Graph;
	FILE *out;
	FILE *report;
	float highest;
	time_t now;
	int used_time;
	int m,ms;
	FILE *mdl_file;
	time_t t;

	B_Options= Options;

/*
#ifdef DNA
   if  ((Block_Init_Sets("Pratt.sets.dna"))==0) {
      printf ("Missing file Pratt.sets.dna\n");
      exit (1);
	}
#else
   if  ((Block_Init_Sets("Pratt.sets"))==0) {
      printf ("Missing file Pratt.sets\n");
      exit (1);
   }
#endif
*/

	if (Options->Input_Symbol_File) {
		if (!Block_Init_Sets(Options->Symbol_File)) {
			fprintf(stderr,"Could not use symbol file %s\n",Options->Symbol_File);
			exit(1);
		}
	}
	else if (!Block_Init_Sets("")) {
		fprintf(stderr,"Error initialising symbol table\n");
		exit(1);
	}

	if (Options->Tree_Input) {
		if ((Tree= Tree_From_File(Options->Tree_Filename))==NULL) {
			fprintf(stderr,"\nERROR READING TREE FILE\n");
			fprintf(stderr,"-----------------------\n");
			fprintf(stderr,"Tree could not be read from file %s\n",Options->Tree_Filename);
			fprintf(stderr,"Please check file and file format\n");
			exit(1);
		}
		Tree_Print (Tree);
		Make_Edge_Sets_Of_Nodes(Tree);
		if (Tree_Make_Seq_Pointers(Tree, seqs, nrseqs)==0)
		{
			printf ("names in tree file does not match sequence names\n");
			exit(1);
		}
	}
	if (Options->Dist_Input)
	{
		int nr_dist;
		nr_dist= read_mat(Dist_Matrix,Options->Dist_Filename);
		if (nr_dist!=nrseqs) {
			printf ("should be the same number of sequences in distance file\n");
			printf ("nr_dist=%d, nrseqs=%d\n",nr_dist,nrseqs);
			exit(1);
		}
	}


/*
	if (Options->Tree_Input && Options->Dist_Input)
		Make_Random_Sequence_Sets(seqs,nrseqs,Tree);

   exit(1);
*/


	Nr_Sets_Block= min(Nr_Sets,Options->Nr_Symbol_Block);

	block= Block_Init(nrseqs,seqs);

	if (Options->Query_Mode) {
		al= Al_From_Sequences(&query_sequence,1);
#ifdef DEBUG
      Alignment_Print(al);
#endif
      Pat_Graph= Pat_Graph_From_Alignment(al,block);
	}
	else if (Options->Alignment_flag) {
      al= Alignment_From_File(Options->Al_Filename);
#ifdef DEBUG
      Alignment_Print(al);
		Alignment_Print_Seqs_Fasta(al,"ut.fasta");
#endif
      Pat_Graph= Pat_Graph_From_Alignment(al,block);
   }
	else if (Options->Use_Short_Sequence_To_Guide_Search) {
		al= Alignment_From_Shortest_Sequence(nrseqs,seqs,Options->Min_Nr_Seqs_Matching);
		assert (al);
#ifdef DEBUG
      Alignment_Print(al);
#endif
      Pat_Graph= Pat_Graph_From_Alignment(al,block);
	}


	printf ("The initial pattern search:\n");

	assert (B_Options->Max_Flex_Prod==Options->Max_Flex_Prod);

	Nr_Sets_In_Search= min(Options->Nr_Symbol_Search1,Nr_Sets_Block);

	Length_Block= block->len;
	if (Allocate_Structures())
		printf ("  -- allocated memory.\n");

	Seqs_Prepare_For_Output(seqs,nrseqs);

	if (B_Options->MDL_flag) {
		if ((mdl_file= fopen("mdl.special","w"))==NULL) {
			fprintf (stderr,"Cannot write to MDL outputfile mdl.special\n");
			exit(1);
		}
	}

	start_search_time= time(NULL);

	if (Hits1= Find_Conserved_Patterns_Using_Graph(block,Pat_Graph)) {
		if (Options->Diagnostic)
			Nr_Entries_Swiss_Prot= Find_Hits_In_SWISS_PROT(Hits1,B_Options->Swiss_Flat_File);

		Hit_List_Find_Statistics(block,Hits1);

		if (Options->Refinement)
		{
			Hit_Print_Hits(block,Hits1,Options->Filename,"a",0);
			Hits2= Refine_Hits(block,Hits1);
		}
		else {
			Hits2= Hits1;
		}

/*
		Hits_Remove_Gappy_Columns (Hits2,block);
*/

		if (Max_Comp_Used(Hits2))  {
			if ((out=fopen(Options->Filename,"a"))!=NULL) {
				fprintf(out,"\n\nWarning: upper limit on pattern length (components) reached\n");
				fprintf(out,"         -- to check if longer patterns may be conserved.\n");
				fprintf(out,"         -- rerun Pratt with increased C and L values\n");
				fclose(out);
			}
			else
				assert(0);
		}
			
  		Hit_Print_Hits(block,Hits2,Options->Filename,"a",1); 
  		Hit_Print_Hits(block,Hits2,Options->Filename,"a",2); 
		if (B_Options->MDL_flag)
        	Hit_Print_mdl_Hits(block,Hits2,mdl_file);

		if ((Hits2->First && (Hits2->First->Pat!=NULL)) && (Options->Print_Seqs_With_Motifs))
		{
			if ((out= fopen(Options->Filename,"a"))==NULL) {
				fprintf (stderr,"Cannot open file %s for writing\n",Options->Filename);
			}
			else {
				fprintf (out,"\n\nPATTERN MATCHES:\n");
				fprintf (out,"each . represents %2d sequence symbols\n",Options->Print_Ratio);
				fprintf (out,"A symbol A-Z,a-z (for example A) in the place of a dot indicates the\n");
				fprintf (out,"starting point of a match to this pattern (in the example; pattern A).\n");
				fprintf (out,"\n");
				Seqs_Print_Out(&out,seqs,nrseqs);
				fclose(out);
			}
		}

	}

	if (Hits2 && Hits2->First && Hits2->First->Pat)
	{
		highest= ((t_Pat_Info*)(Hits2->First->Pat))->information;
		m= ((t_Pat_Info*)(Hits2->First->Pat))->nr_match;
		ms= ((t_Pat_Info*)(Hits2->First->Pat))->nr_match_swiss;
	}
	else
	{
		highest= 0.0;
		m=ms=0;
	}
	if (B_Options->MDL_flag)
		fclose(mdl_file);

	now= time(NULL);
	/*
	used=difftime(now,start_sec);
	*/
	used_time= (now-start_sec);

	printf ("Total running time: %4d seconds%s\n",used_time,(Truncated?" (Truncated)":""));

	report= fopen("report","a");
	if (Hits2&& Hits2->First && Hits2->First->Pat) {
		fprintf (report,"%40s %4d %8d %8d %8.4f ", Options->Filename,used_time,sum_length_sequences,Graph_Nr_Nodes,((t_Pat_Info*)(Hits2->First->Pat))->information);
		Pat_Print(report,((t_Pat_Info*)(Hits2->First->Pat))); 
		fprintf (report,"\n");
	}
	else 
		fprintf (report,"%40s %4d %8d %8d    0.0\n", Options->Filename,used_time,sum_length_sequences,Graph_Nr_Nodes);
	fclose(report);

/*
	if (Hits2&& Hits2->First && Hits2->First->Pat) {
		fprintf (report,"%40s %8.4f %6.4f %6.4f  ", Options->Filename,((t_Pat_Info*)(Hits2->First->Pat))->information, ((t_Pat_Info*)(Hits2->First->Pat))->sensitivity, ((t_Pat_Info*)(Hits2->First->Pat))->PPV);
		Pat_Print(report,((t_Pat_Info*)(Hits2->First->Pat))); 
		fprintf (report,"\n");
	}
	else 
		fprintf (report,"%40s\n ", Options->Filename);
*/


/*
	code used for the blocks project:
*/
/*
	if (Hits2&& Hits2->First && Hits2->First->Pat) {
		fprintf (report,"%40s %d %8.4f %d %d ", Options->Filename,(int) used, ((t_Pat_Info*)(Hits2->First->Pat))->information, nrseqs, ((t_Pat_Info*)(Hits2->First->Pat))->nr_match);
		Pat_Print(report,((t_Pat_Info*)(Hits2->First->Pat))); 
		fprintf (report,"\n");
	}
	else 
		fprintf (report,"%40s\n ", Options->Filename);

	fclose(report);
*/


	Seqs_Prepare_For_Output(seqs,nrseqs);

	if (Hits2->First && (Hits2->First->Pat!=NULL)) {
		if ((out= fopen(Options->Filename,"a"))==NULL) {
			printf ("Pratt cannot open file %s for writing \n",Options->Filename);
			printf ("Change write access or use another filename and run Pratt again\n");
			printf ("Ciao\n");
		}
		fprintf (out,"Number of patterns evaluated by Pratt:%d\n",Nr_Patterns_Evaluated);
		fprintf (out,"Total running time: %4d seconds %s\n",used_time,(Truncated?"(Truncated)":""));
		fprintf (out,"\n\n");
		fclose(out);

		printf ("number of patterns evaluated:%d\n",Nr_Patterns_Evaluated);
		printf ("\n");
		printf ("Results printed to file %s\n",B_Options->Filename);
		printf ("\n");
	} 
	else {
		if ((out= fopen(Options->Filename,"a"))==NULL) {
			printf ("Pratt cannot open file %s for writing \n",Options->Filename);
			printf ("Change write access or use another filename and run Pratt again\n");
			printf ("Ciao\n");
		}
		fprintf(out,"No patterns found\n");
		fclose(out);
		printf("\nNo patterns found\n");
	}
}



void flush_header()
/******************************************************************/
{
	printf ("\nPattern search statistics:\n\n");
	printf ("  Proportion  | Number of | Score of best \n");
	printf ("   of search  |  patterns | pattern found \n");
	printf (" (estimated)  | looked at | so far    \n");
	printf (" -------------+-----------+--------------\n");
}


void flush_status(t_Hits* Hits)
/******************************************************************/
{
/*
	time_t now;
	double used;
	int left;
	double f;

	now= time(NULL);
	used=difftime(now,start_search_time);

	f= (double)nr_now/(double)nr_total;
	left= (int) ((1.0-f)/(f/used));

	if (Hits->First->Information<0.1)
		printf ("\r %9.0f%%   |%10d | none    | %4d secs",100.0*(float)nr_now/(float)nr_total,Nr_Patterns_Evaluated,left);
	else
		printf ("\r %9.0f%%   |%10d | %10.4f  | %4d secs",100.0*(float)nr_now/(float)nr_total,Nr_Patterns_Evaluated,Hits->First->Information,left);
*/
	if (Hits->First->Information<0.1)
		printf ("\r %9.0f%%   |%10d | none    ",100.0*(float)nr_now/(float)nr_total,Nr_Patterns_Evaluated);
	else
		printf ("\r %9.0f%%   |%10d | %10.4f  ",100.0*(float)nr_now/(float)nr_total,Nr_Patterns_Evaluated,Hits->First->Information);

	fflush(stdout);
}


static int Calculate_Cardinality(t_Pattern *Pat, int refinement_flag)
/******************************************************************/
{
	static int i;
	static unsigned int *vec1;
	static unsigned int *end1;
	static unsigned int *vec2;
	static unsigned int* temp_vector=NULL;
	static size_t length_temp_vector_bytes;
	static int temp;
	static int num_sub;

	if (refinement_flag)
		num_sub= Pat->num_sub_new;
	else
		num_sub= Pat->num_sub;

	if (temp_vector==NULL) {
		length_temp_vector_bytes=sizeof(unsigned int)*Length_Block;
		temp_vector= malloc(length_temp_vector_bytes);
	}

	temp= 0;
	if (num_sub==1)
	{
		vec1= Pat->bit_vecs[0];
		end1= vec1+Length_Block;
		while (vec1<end1) {
			temp+=NUM_ONES32(*vec1);
			vec1++; 
		}
	}
	else
	{
		memcpy(temp_vector,Pat->bit_vecs[0],length_temp_vector_bytes);
		for (i=1;i<num_sub;i++) {
			vec1= Pat->bit_vecs[i];
			end1= vec1+Length_Block;
			vec2= temp_vector;
			while (vec1<end1)
				*vec2++|= *vec1++;
		}
		vec1= temp_vector;
		end1= temp_vector+Length_Block;
		while (vec1<end1) {
			temp+=NUM_ONES32(*vec1);
			vec1++; 
		}
	}
	return temp;
}

/*
	This file is part of the Pratt program source code
	Copyright 1996 Inge Jonassen, Dept. of Informatics
	University of Bergen. 
	email: inge@ii.uib.no
	More information on Pratt: 
	http://www.ii.uib.no/~inge/Pratt.html
*/
