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
#include "sequence.h"
#include "hit.h"
#include "al.h"
#include "pattern.h"
#define BLOCK_MAIN
#include "block.h"
#include "menu.h"
#include "scan.h"
#include "swiss.h"
#include "tree.h"
#include "mst.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))


 
/*
#define DEBUG
*/


/**************** VARIABLES GLOBAL TO THIS FILE *************/
static float Freq[26];
static int Swiss_matched= 0;
static unsigned int *temp_vector;
static unsigned int *temp_vector2;


/**************** GLOBAL VARIABLES DEFINED HERE *************/
unsigned int Member_In_Set[MAX_NR_SETS];
int Number_In_Set[MAX_NR_SETS];
float Freq_Set[MAX_NR_SETS];
float Bit_Set[MAX_NR_SETS];
int Nr_Sets;
int Nr_Sets_Block;
int Nr_Sets_Unit;
int Nr_Sets_In_Search;
char AA_Tab[20];
int AA_Index[26];
int Num_Ones[LOOK_UP_ONES];
int Length_AA_Set;
unsigned int AA_In_Set[26][10];
t_Generalize_Symbol_List **Gen_Symbol_Lists;
int sum_length_sequences;

float apr_single[26]= 
	{ .08713,  /* A */ 
	  .0    ,  /* B */
	  .03347,  /* C */
	  .04687,  /* D */
	  .04953,  /* E */
	  .03977,  /* F */
	  .08861,  /* G */
	  .03362,  /* H */
	  .03689,  /* I */
	  .0,      /* J */
	  .08048,  /* K */
	  .08536,  /* L */
	  .01475,  /* M */
	  .04043,  /* N */
	  .0,      /* O */ 
	  .05068,  /* P */
	  .03826,  /* Q */
	  .04090,  /* R */
	  .06958,  /* S */
	  .05854,  /* T */
	  .0,      /* U */
	  .06472,  /* V */
	  .01049,  /* W */
	  .0,      /* X */
	  .02992,  /* Y */
	  .0};     /* Z */

float Shannon;
/**************** GLOBAL VARIABLES DEFINED ELSEWHERE ********/
extern t_sequence** seqs;
extern int nrseqs;
extern int virt_nrseqs;
extern t_Block_Options *B_Options;
extern int Nr_Entries_Swiss_Prot;
extern t_tree *Tree;
extern float Dist_Matrix[MAX_NUM_DIST][MAX_NUM_DIST];

/**************** LOCAL ROUTINES ****************************/
static int Block_Init_Generalize_List();


static char PrattSets[200][30];

static char PrattSetsDefault[][30]=
{"I",
"L",
"V",
"M",
"F",
"Y",
"W",
"H",
"K",
"T",
"C",
"E",
"A",
"P",
"G",
"S",
"D",
"N",
"Q",
"R",
"RKH",
"DE",
"NQ",
"ST",
"ILV",
"FWY",
"AG",
"DERKH",
"DERK",
"TSND",
"SND",
"AGS",
"HWYF",
"SNDEQ",
"SNDEQR",
"PTSNDEQ",
"TSNDEQ",
"PTSND",
"AGTSND",
"PAGTSND",
"AGTSNDEQ",
"VCAGTSND",
"VCAGTS",
"VCAGT",
"VCAG",
"VCAGP",
"VCAGPT",
"LIVAGP",
"FMLIV",
"\0"};


/*
static char PrattSetsDefault[][30]=
{"I",
"L",
"V",
"M",
"F",
"Y",
"W",
"H",
"K",
"T",
"C",
"E",
"A",
"P",
"G",
"S",
"D",
"N",
"Q",
"R",
"RKH",
"DE",
"NQ",
"ST",
"ILV",
"FWY",
"AG",
"DERKH",
"DERK",
"TSND",
"SND",
"AGS",
"HWYF",
"SNDEQ",
"SNDEQR",
"SNDEQRKH",
"SNDEQRKHP",
"PTSNDEQRKH",
"TSNDEQRKH",
"PTSNDEQ",
"TSNDEQ",
"PTSND",
"AGTSND",
"PAGTSND",
"AGTSNDEQ",
"PAGTSNDEQ",
"AGTSNDEQRK",
"PVCAGTSND",
"VCAGTSND",
"VCAGTS",
"VCAGT",
"VCAG",
"VCAGP",
"VCAGPT",
"LIVCAGPT",
"LIVAGP",
"LIV",
"FMLIV",
"FMLIVCAG",
"PFMLIVCAG",
"FMLIVCAGT",
"PFMLIVCAGT",
"FMLIVCAGTK",
"HWYFMLIV",
"HWYFM",
"RKHWYF",
"QERKHWY",
"QERKHYFM",
"LIVMFYWC"
"\0"};

*/




#ifdef DEBUG
#define CHECK(a) assert(a)
#else
#define CHECK(a)
#endif



static int Block_Init_Sets_From_Array(char sets[][30])
/******************************************************************/
{
	char Line[100];
	int h;
	int i,j;
	int unit_ended;

/*
	printf ("  Initialising internal block data structure for pattern matching\n");
*/
	
	for (i=0;i<MAX_NR_SETS;i++)
		Number_In_Set[i]= 0;

	unit_ended= 0;
	Nr_Sets= Nr_Sets_Unit= 0;
	h=0;
	while (sets[h][0]!='\0')
	{
		char *c;
		char letter;
		strcpy (Line,sets[h]);
		assert (Nr_Sets<MAX_NR_SETS);
		Member_In_Set[Nr_Sets]= 0;
		for (c= Line; (*c); c++)
		{
			if (isupper(*c)) {
				Member_In_Set[Nr_Sets]|= 1<< ((*c)-'A');
				letter= (*c);
				Number_In_Set[Nr_Sets]++;
			}
		}
		if (Number_In_Set[Nr_Sets]==1)
		{
			AA_Tab[Nr_Sets_Unit]= letter;
			AA_Index[letter-'A']= Nr_Sets_Unit;
			Nr_Sets_Unit++;
			if (unit_ended) {
			   printf ("Wrong order in Pratt.sets\n");
			   printf ("The single amino acid lines should come first\n");
			   exit (1);
 			 }
		}

		if (Number_In_Set[Nr_Sets]<20) {
			if (Number_In_Set[Nr_Sets]>1)
				unit_ended= 1;
			Nr_Sets++;
		}
		else {
			Number_In_Set[Nr_Sets]= 0;
			Member_In_Set[Nr_Sets]= 0;
		}	
		h++;
	}



	if (Nr_Sets<B_Options->Nr_Symbol_Block) {
		B_Options->Nr_Symbol_Block=B_Options->Nr_Symbol_Search1=Nr_Sets;
	}

	printf ("  %d of the %d potential pattern symbols are used in the initial pattern search\n",B_Options->Nr_Symbol_Block,Nr_Sets);

/*
	for (i=0;i<LOOK_UP_ONES;i++) {
		int num= 0;
		int k;
		for (k=0;k<16;k++) {
			if (i&(1<<k))
				num++;
		}
		Num_Ones[i]= num;
	}
*/

/* John Collins' code: */
{
  int *Num_Ptr1, *Num_Ptr2, *Num_Ptr_End;
  int one = 1;
  int POTWO = 1;
  
  Num_Ones[0] = 0;
  while( POTWO < LOOK_UP_ONES) {
    Num_Ptr1 = Num_Ones + POTWO;
    *Num_Ptr1++ = one;
    Num_Ptr2 = Num_Ones + one;
    Num_Ptr_End = Num_Ptr1 + POTWO;
    while( Num_Ptr1 < Num_Ptr_End) {
      *Num_Ptr1++ = *Num_Ptr2++ + one ;
      }
    POTWO <<= one;
    }
}
/*
	for (i=0;i<LOOK_UP_ONES;i++) {
		assert (Num_Ones[i]==Num_Ones2[i]);
	}
*/
			



	Length_AA_Set= Nr_Sets/PER_WORD + 1;

	for (i=0;i<26;i++)
		for (j=0;j<Length_AA_Set;j++)
			AA_In_Set[i][j]= 0;
	
	for (i=0;i<26;i++)
	{
		for (j=0;j<Nr_Sets;j++)
		{
			if (Member_In_Set[j]&(1<<i))
				AA_In_Set[i][j/PER_WORD]|= (1<<(j%PER_WORD));
		}
	}

	Block_Init_Generalize_List();

	return 1;
}


int Block_Init_Sets(char *filename)
/******************************************************************/
{
	FILE* file;
	char Line[100];
	int i,j;
	int unit_ended;

	printf ("Initialising the set of possible pattern symbols and the block data structure\n");

	if (strlen(filename)==0) {
		printf ("  Using default pattern symbols:\n");
		printf ("    if alternative pattern symbols are wanted,\n");
		printf ("    please compile your own Pratt.sets file.\n");
		return Block_Init_Sets_From_Array(PrattSetsDefault);
	}

	if ((file= fopen(filename,"r"))==NULL) {
		fprintf(stderr,"Could not open file %s for reading\n",filename);
		return 0;
	}

	printf ("  Using symbols in file %s\n",filename);
	i=0; 
	while ((fgets(Line,100,file))!=NULL) {
		strcpy (PrattSets[i],Line);
		i++;
	}
	sprintf(PrattSets[i],"\0");
	fclose(file);
	return Block_Init_Sets_From_Array(PrattSets);
}


static int Block_Init_Generalize_List()
/******************************************************************/
{
	int i,j,k;
	unsigned int symbol_set;
	t_Generalize_Symbol_List* list;

	Gen_Symbol_Lists= malloc(B_Options->Nr_Symbol_Block*sizeof(t_Generalize_Symbol_List*));

	for (i=0;i<B_Options->Nr_Symbol_Block;i++)
	{
		list= Gen_Symbol_Lists[i]= malloc(sizeof(t_Generalize_Symbol_List));
		list->nr_symbols= 0;
		symbol_set= Member_In_Set[i];
		for (j=0;j<B_Options->Nr_Symbol_Block;j++)
		{
			if ((Member_In_Set[j]&symbol_set)==symbol_set)
			{
				list->symbols[list->nr_symbols++]= j;
				if (list->nr_symbols>=MAX_NR_GEN)
				{
					printf ("\n*******************\nPratt fails:\n*****************\n");
					printf ("Too many ambiguous symbols allowed during initial search\n");
					printf ("gives too many possible generalisations of the pattern symbol ");
					for (k=0;k<30;k++)
						if (symbol_set & (1<<k))
							printf ("%c",(char)(k+'A'));
					printf ("\n");	
					printf ("current limit is %d.\n",MAX_NR_GEN);
					printf ("Reduce the amount of ambiguity to be allowed or increase MAX_NR_GEN\n");
					printf ("in block.h and recompile everyting.\n");
				}
				assert (list->nr_symbols<MAX_NR_GEN);
			}
		}
	}
}


t_block* Block_Init(int nr_seqs, t_sequence *Seqs[])
/******************************************************************/
{
	t_block* block;
	int vector_len;
	int i,j,k,kk;
	int len;
	unsigned int *vec1;
	float *fpt1,*fpt2;
	int average_length;
	int restricted;
	unsigned int *allowed_vector;;

	printf ("  Filling in the block data structure\n");

	block= malloc(sizeof(t_block));
	assert (block);

/*
	if (nr_seqs==1) {
		vector_len= Seqs[0]->length/PER_WORD + ((Seqs[0]->length)%PER_WORD ? 1 : 0 );
		block->len= vector_len;
		block->start=  malloc ((vector_len+1)*sizeof(int));
		block->length=  malloc ((vector_len+1)*sizeof(int));
		block->seq=  malloc ((vector_len+1)*sizeof(int));
		assert (block->start && block->length && block->seq);

		for (i=0;i<vector_len;i++) {
			block->seq[i]= i;
			block->start[i]=i;
			block->length[i]=1;
		}
		block->start[vector_len]=vector_len;
		virt_nrseqs= block->len;
	}
	else 
*/
	{
		virt_nrseqs= nr_seqs;

		block->start=  malloc ((nr_seqs+1)*sizeof(int));
		block->length=  malloc (nr_seqs*sizeof(int));
		assert (block->start && block->length);
	
		if (!block || !(block->start) || !(block->length))
			return NULL;
	
		average_length= 0;
		len= 0;
		block->start[0]= 0;
		for (i=0;i<nr_seqs;i++)
		{
			average_length+= Seqs[i]->length;
			len+= Seqs[i]->length;
			block->length[i]= Seqs[i]->length/PER_WORD + 
                        	((Seqs[i]->length)%PER_WORD ? 1 : 0 );
			block->start[i+1]= block->start[i]+block->length[i];
		}
	
		sum_length_sequences= average_length;
	
		average_length/= nr_seqs;
	
		printf ("    %d sequences of average length %d\n",nr_seqs,average_length);
	
		vector_len= block->start[nr_seqs];
	
		block->len= vector_len;
	
		block->seq= (int*) malloc (block->len*sizeof(int));
		assert (block->seq);
	
		if (!(block->seq))
			return NULL;
	
		for (i=0;i<nr_seqs;i++)
			for (j=block->start[i];j<block->start[i]+block->length[i];j++)
				block->seq[j]= i;
	}

	block->vectors= malloc(B_Options->Max_Length*Nr_Sets_Block*sizeof(unsigned int*));
	assert (block->vectors);

	for (i=0; i<B_Options->Max_Length*Nr_Sets_Block;i++) {
		block->vectors[i]= calloc(vector_len,sizeof(unsigned int));
		assert (block->vectors[i]);
	}

	if (B_Options->Restrictions_Input==complete) {
		restricted= 1;
		allowed_vector= calloc(block->len,sizeof(unsigned int));
		assert (allowed_vector);
		for (i=0;i<nr_seqs;i++) {
			if (Seqs[i]->restricted) {
				for (j=0;j<Seqs[i]->restricted;j++)
					for (k=Seqs[i]->first[j];k<=Seqs[i]->last[j];k++)
						allowed_vector[block->start[i]+(k/PER_WORD)]|= 1<<(k%PER_WORD);
			}
			else {
			  for (k=0;k<Seqs[i]->length;k++)
				  allowed_vector[block->start[i]+(k/PER_WORD)]|= 1<<(k%PER_WORD);
		  }
		}
	}
	else
		restricted= 0;


	for (i=0;i<nr_seqs;i++) {
		char *cpt1,*cpt2; 
		int seq_length;
		cpt1= Seqs[i]->seq;
		cpt2= cpt1+Seqs[i]->length;
		seq_length= Seqs[i]->length;
		for (j=0;j<seq_length;j++) {
			int a;
			int l;

			a= (*cpt1++)-'A';
			CHECK (a>=0 && a<26);
			if (a<0 || a>25)
				continue;
			/* if complete restrictions then position j in sequence i has to be
				within allowed window */
			if (restricted &&
				 (Seqs[i]->restricted>0) &&
				 ((allowed_vector[block->start[i]+(j/PER_WORD)] & (1<<(j%PER_WORD)))==0) )
				continue;
			for (l=0;l<Nr_Sets_Block;l++) {
				if ((1<<a) & Member_In_Set[l]) 
				{
					for (k=0;k<=min(j,B_Options->Max_Length-1);k++) {
						for (kk= 0;kk<=B_Options->Sloppy;kk++) {
							if ((k+kk)<0)
								continue;
							if ((k+kk)>=B_Options->Max_Length-1)
								continue;
							if (j-(k+kk)<0)
								continue;
						/* Position k(+kk) within block j-k in sequence i is in set l */
						if (nr_seqs>1)
							block->vectors[Nr_Sets_Block*(k)+l][block->start[i]+((j-(k+kk))/PER_WORD)]|=
									1<<((j-(k+kk))%PER_WORD);
						else
							block->vectors[Nr_Sets_Block*(k)+l][((j-(k+kk))/PER_WORD)]|= 1<<((j-(k+kk))%PER_WORD);
						}
					}
				}
			}
		}
	}

	fpt1= Freq; fpt2= Freq+26;
	while (fpt1<fpt2)
		*fpt1++= 0.0;

	for (i=0;i<nr_seqs;i++)
	{
		char *cpt1,*cpt2; 
		cpt1= Seqs[i]->seq;
		cpt2= cpt1+Seqs[i]->length;
		while (cpt1<cpt2)
			Freq[(*cpt1++)-'A']+= 1.0;
	}

	fpt1= Freq; fpt2= Freq+26;
	while (fpt1<fpt2)
		*fpt1++/= (float)len;

	fpt1= apr_single; fpt2= apr_single+26;
	Shannon= 0.0;
	while (fpt1<fpt2) {
		if (*fpt1>0.0)
			Shannon-= *fpt1*log(*fpt1)/log(2.0);
		fpt1++;
	}

	for (i=0;i<Nr_Sets;i++) {
		Freq_Set[i]= 0.0;
		for (j=0;j<26;j++) {
			if (Member_In_Set[i]&(1<<j))
				Freq_Set[i]+= apr_single[j];
		}
		Bit_Set[i]= Shannon;
		for (j=0;j<26;j++) {
			if (Member_In_Set[i]&(1<<j))
			{
				float fj= apr_single[j]/Freq_Set[i];
				Bit_Set[i]+= fj*log(fj)/log(2.0);
			}
		}
	}

/*
	for (i=0;i<Nr_Sets_Unit;i++) 
	{
		printf ("%f,%f\n", Freq_Set[i],apr_single[AA_Tab[i]-'A']);
		assert(Freq_Set[i]==apr_single[AA_Tab[i]-'A']);
	}
*/

	temp_vector= malloc(block->len*sizeof(unsigned int));
	temp_vector2= malloc(block->len*sizeof(unsigned int));
	assert (temp_vector && temp_vector2);

	printf ("-- finished initialising internal block data structure.\n\n");

	return block;
}

void Block_Print_Hits_To_Sequence(char symbol,t_block* block, t_Pat_Info *Pat_Info)
/******************************************************************/
{
	unsigned int* stat;
	int length;
	int i,j,k;


	stat= Pat_Info->bit_vector;
	length= Pat_Info->length;

	for (i=0;i<nrseqs;i++) {
		for (j=0;j<block->length[i];j++) {
			for (k=0;k<PER_WORD;k++) {
				if (stat[block->start[i]+j] & (1<<k) ) {
					Seq_Add_Symbol_To_Output(seqs[i],(int)(j*PER_WORD+k),symbol);
				}
			}
		}
	}

}


void Block_Print_Hits(FILE* file,t_block* block, t_Pat_Info *Pat_Info)
/******************************************************************/
{
	unsigned int* stat;
	int length;
	int i,j,k;

	Special_Print_Pattern_Matches(file,block,seqs,nrseqs,Pat_Info);
}

int Pat_Info_Divergence(t_Pat_Info *Pat,t_block *block)
/******************************************************************/
{
	int Inc[MAX_NUM_DIST];
	int i,j,k;
	int found;
	unsigned int* bit_vector;

	if (B_Options->Tree_Input || B_Options->Dist_Input) 
	{
		bit_vector= Pat->bit_vector;

		for (i=0;i<nrseqs;i++) {
			found= 0;
			for (j=0;j<block->length[i];j++) {
				for (k=0;k<PER_WORD;k++) {
					if (bit_vector[block->start[i]+j] & (1<<k) ) {
						found=1;
					}
				}
			}
			Inc[i]= found;
		}
	}

	if (B_Options->Tree_Input) {
		Pat->seq_hit_divergence= Divergence(Tree,Inc,seqs,nrseqs);
		Pat->fitness= Pat->information*Pat->seq_hit_divergence;
		return 1;
	}
	else {
		Pat->seq_hit_divergence= 0.0;
	}	

	if (B_Options->Dist_Input) {
		char mst_tree[1000];
		Pat->seq_hit_div_mst= mst_weight(Dist_Matrix,nrseqs,Inc,mst_tree);
		Pat->fitness= Pat->information*Pat->seq_hit_div_mst;
		return 1;
	}
	else {
		Pat->seq_hit_div_mst= 0.0;
	}

	return 0;
}


void Statistics_Refined_Pattern(t_block* block, t_Pat_Info *Pat)
/******************************************************************/
{
	unsigned int* bit_vector;
	int i,j,k;
	int num, num_seq;
	int found;
	float TP,TN,FP,FN;


	if (Pat==NULL) {
		printf ("Statistics_Refined_Pattern NULL Pat_Info\n");
		exit (1);
		return;
	}

	bit_vector= Pat->bit_vector;

	num= num_seq= 0;

	for (i=0;i<virt_nrseqs;i++) {
		found= 0;
		for (j=0;j<block->length[i];j++) {
			for (k=0;k<PER_WORD;k++) {
				if (bit_vector[block->start[i]+j] & (1<<k) ) {
					found=1;
					num++;
				}
			}
		}
		if (found)
			num_seq++;
	}
	Pat->nr_occ= num;
	Pat->nr_match= num_seq;
	/*
	assert (num>0 && num_seq>0);
	*/

	if (B_Options->MDL_flag)
		Pat->info_MDL= Pat->information - B_Options->MDL_constant0 - ((B_Options->MDL_constant1*Pat_Info_Num_Char(Pat)+B_Options->MDL_constant2*Pat_Info_Num_X(Pat)+B_Options->MDL_constant3)/((float)max(1,num_seq))); 
	else
		Pat->info_MDL= Pat->information;

	if (Pat_Info_Divergence(Pat,block)) 
		;
	else if (B_Options->Diagnostic)
	{
		TP= num_seq;              /* is in family and matches pattern */
		FN= nrseqs-num_seq;       /* is in family, but dont match p.  */
		FP= Pat->nr_match_swiss-num_seq;  /* matches p, not in family */
		TN= Nr_Entries_Swiss_Prot-(nrseqs+FP);
								  		/* not in family, dont match        */

		Pat->correlation= (TP*TN-FP*FN)/sqrt((TP+FP)*(FP+TN)*(TN+FN)*(FN+TP));
		Pat->sensitivity= TP/(TP+FN);
		Pat->specificity= TN/(TN+FP);
		Pat->fitness= Pat->PPV= TP/(TP+FP);
	}
	else
	{
		Pat->sensitivity= 0.0;
		Pat->specificity= 0.0;
		if (B_Options->MDL_flag)
			Pat->fitness= Pat->info_MDL; /* rank patterns according to MDL measure */
		else
			Pat->fitness= Pat->information;
	}
}



void Hit_List_Find_Statistics(t_block* block, t_Hits* hits)
/******************************************************************/
{
	t_Hit_Entry *Entry;
	t_Pat_Info *Pat;
	unsigned int* bit_vector;
	int i,j,k;
	int num, num_seq;
	int found;
	float TP,TN,FP,FN;

	int Inc[MAX_NUM_DIST];

	for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next)
	{
		if ((Pat=Entry->Pat)==NULL)
			continue;

		bit_vector= Pat->bit_vector;

		num= num_seq= 0;

		for (i=0;i<virt_nrseqs;i++) {
			found= 0;
			for (j=0;j<block->length[i];j++) {
				for (k=0;k<PER_WORD;k++) {
					if (bit_vector[block->start[i]+j] & (1<<k) ) {
						found=1;
						num++;
					}
				}
			}
			if (found) {
				num_seq++;
				Inc[i]= 1;
			}
			else
				Inc[i]= 0;
		}
		Pat->nr_occ= num;
		Pat->nr_match= num_seq;

		if (Pat_Info_Divergence(Pat,block))
			continue;

		if (!B_Options->Diagnostic) {
			Pat->fitness= Pat->information;
			continue;
		}

		TP= num_seq;              /* is in family and matches pattern */
		FN= nrseqs-num_seq;       /* is in family, but dont match p.  */
		FP= Pat->nr_match_swiss-num_seq;  /* matches p, not in family */
		TN= Nr_Entries_Swiss_Prot-(nrseqs+FP);
										  /* not in family, dont match        */

		Pat->correlation= (TP*TN-FP*FN)/sqrt((TP+FP)*(FP+TN)*(TN+FN)*(FN+TP));
		Pat->sensitivity= TP/(TP+FN);
		Pat->specificity= TN/(TN+FP);
		Pat->PPV= TP/(TP+FP);
		Pat->fitness=Pat->PPV;
	}
}

void Hit_Print_mdl_Hits(t_block* block,t_Hits* hits,FILE *mdl_file)
/******************************************************************/
{
	t_Hit_Entry *Entry;
	t_Pat_Info *Pat_Info;
	char symbol;
	int i, found,j,k;
	unsigned int* bit_vector;
	int index;

	assert (B_Options->MDL_flag);

	index=1;
	for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next)
	{
		if (!Entry->Pat)
			continue;

		Pat_Info= Entry->Pat;

		fprintf (mdl_file,"%10.5f %10.5f %d ", Pat_Info->info_MDL,Pat_Info->information,index);
		Pat_Print(mdl_file,Pat_Info);
	
		bit_vector= Pat_Info->bit_vector;
	
		for (i=0;i<nrseqs;i++) {
			found= 0;
			for (j=0;j<block->length[i];j++) {
				for (k=0;k<PER_WORD;k++) {
					if (bit_vector[block->start[i]+j] & (1<<k) ) {
						found=1;
					}
				}
			}
			if (found) 
				fprintf (mdl_file," %d",i);
		}
		fprintf (mdl_file,"\n");
		index++;
	}

/*
	fprintf (mdl_file,"\n");
*/
}

void Hit_Print_Hits(t_block* block,t_Hits* hits,char *filename, char* mode,int hit_flag)
/* hit_flag = 0 : before refinement -
   hit_flag = 1 : after refinement
   hit_flag = 2 : after refinement, print aligments
*/
/******************************************************************/
{
	t_Hit_Entry *Entry;
	FILE* Pattern_File;
	int rank;
	int matches;
	t_Pat_Info *Pat_Info;
	char symbol;
	char s1[100],s2[100],s3[100],s4[100];

	Pattern_File= fopen(filename,mode);
	
	if (Pattern_File==NULL)
	{
		printf ("Hit_Print_Hits %s %d: ", __FILE__,__LINE__);
		printf ("Could not write to file %s\n",filename);
		return;
	}

	if (hits->First->Pat==NULL)
		return;

	if (hit_flag==0)
		printf ("\n\nBest Patterns before refinement:\n");
	else if (hit_flag==1)
		printf ("\n\nBest Patterns %s:\n",((B_Options->Refinement==1) ? "(after refinement phase)" : ""));
	else 
		printf ("Occurrences in sequences are being written to file %s\n\n",filename);

	if (hit_flag==0)
		fprintf (Pattern_File,"\n\nBest Patterns before refinement:\n");
	else if (hit_flag==1)
		fprintf (Pattern_File,"\n\nBest Patterns %s:\n",((B_Options->Refinement==1) ? "(after refinement phase)" : ""));
	else 
		fprintf (Pattern_File,"\n\nBest patterns with alignments:\n");

/*
	fprintf (Pattern_File,"         corr sens spec   fitness       div      mst        family      swiss\n");
*/
	fprintf (Pattern_File,"       %s  fitness   %s %s hits(seqs)  %s  Pattern\n",
		(B_Options->Diagnostic?"  corr sens spec":""),
	   (B_Options->Tree_Input?"  div. ":""),
		(B_Options->Dist_Input?"  mcst ":""),        
		(B_Options->Diagnostic?"    SWISS ":""));

	if (hit_flag!=2) {
		printf ("       %s  fitness   %s %s hits(seqs)  %s  Pattern\n",
			(B_Options->Diagnostic?"  corr sens spec":""),
	   	(B_Options->Tree_Input?"  div. ":""),
			(B_Options->Dist_Input?"  mcst ":""),        
			(B_Options->Diagnostic?"    SWISS ":""));
	}

	symbol='A'; rank= 1; 

	/* do not print symbols beofre refinement if refinement is done */
	if (hit_flag==0) 
		symbol= 0;

	for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next)
	{
		if (!Entry->Pat)
			continue;

		Pat_Info= Entry->Pat;

		if (hit_flag==2 && rank>B_Options->Number_Of_Occ_Output)
			/* do not print more than B_Options->Number_Of_Occ_Output alignments */
			;
		else {
			if (B_Options->Diagnostic){
			/*
				sprintf (s1,"%6.4f %4.2f %4.2f",Pat_Info->correlation,Pat_Info->sensitivity,Pat_Info->specificity);
			*/
				sprintf (s1,"%6.4f %4.2f %4.2f",Pat_Info->correlation,Pat_Info->sensitivity,Pat_Info->PPV);
				sprintf (s4,"%4d(%4d)",Pat_Info->nr_occ_swiss, Pat_Info->nr_match_swiss);
			}
				else {
				sprintf(s1,"");
				sprintf(s4,"");
			}

			if (B_Options->Tree_Input)
				sprintf(s2,"%8.4f",Pat_Info->seq_hit_divergence);
			else	
				sprintf(s2,"");

			if (B_Options->Dist_Input)
				sprintf(s3,"%8.4f", Pat_Info->seq_hit_div_mst);
			else	
				sprintf(s3,"");

			fprintf (Pattern_File,"%c %3d: %s %8.4f %s %s %4d(%4d) %s  ",
	  	 		(symbol!=0?symbol:' '),rank,s1, Pat_Info->fitness,s2,s3,Pat_Info->nr_occ, Pat_Info->nr_match,s4);

			Pat_Print(Pattern_File,Pat_Info); fprintf (Pattern_File,"\n");

			if (hit_flag!=2) {
				printf ("%c %3d: %s %8.4f %s %s %4d(%4d) %s  ",
	  	 		(symbol!=0?symbol:' '),rank,s1, Pat_Info->fitness,s2,s3,Pat_Info->nr_occ, Pat_Info->nr_match,s4);

			Pat_Print(stdout,Pat_Info); printf ("\n");
			}
		}
	
		if (hit_flag==2) {
			if (symbol!=0)
				Block_Print_Hits_To_Sequence(symbol,block,Pat_Info);

			if (rank<=B_Options->Number_Of_Occ_Output) {
				Block_Print_Hits(Pattern_File,block,Pat_Info);
				fprintf (Pattern_File,"\n");
			}
		}

		rank++;

		if (symbol=='Z')
			symbol='a';
		else if (symbol=='z')
			symbol= 0;
		else if (symbol!=0)
			symbol++;
	}
	printf ("\n");

	fclose (Pattern_File);
}

static unsigned int **Seq_Hit_Vector=NULL;
static int Length_Seq_Hit_Vector;

void Matching_Prepare(int nrseqs, int max_flex)
/******************************************************************/
{
	int i;

	Length_Seq_Hit_Vector = 1+(nrseqs/PER_WORD);

	Seq_Hit_Vector= malloc((2+max_flex*2)*sizeof(unsigned int*));
	assert (Seq_Hit_Vector);
	for (i=0;i<2+2*max_flex;i++) {
		Seq_Hit_Vector[i]= calloc(Length_Seq_Hit_Vector,sizeof(unsigned int));
		assert (Seq_Hit_Vector[i]);
	}
}


unsigned int Nr_Match(int s1, int s2, int min_dist, int max_dist, int max_flex, t_block *block)
/******************************************************************/
/* here min_dist, max_dist, min_d, max_d is difference in position*/
/* that is c-x(0,1)-c -> min_dist=min_d=1, max_dist=2, max_d=3    */
/******************************************************************/
{
	int i,j,k,l,m,n;
	int last;
	int seq_index;
	int nr;
	int min_d, max_d;
	unsigned int mask;
	unsigned int seq_mask;
	int act_d;
	int flex;
	int num_seq;

	max_dist++; min_dist++;

	assert (Seq_Hit_Vector);

	min_d= max(1,max_dist-max_flex);
	max_d= min(min_dist+max_flex,B_Options->Max_Length-1);

	flex= max_dist-min_dist;

	for (i=min_d;i<=max_d;i++)
	{
		memset(Seq_Hit_Vector[i-min_d],0,Length_Seq_Hit_Vector*sizeof(unsigned int));
		j= 0;
		while (j<block->len) 
		{
			if ((block->vectors[s1][j] & block->vectors[Nr_Sets_Block*i+s2][j])!=0)
			{
				seq_index= block->seq[j];
				assert (i-min_d>=0 && i-min_d <=2*max_flex);
				Seq_Hit_Vector[i-min_d][seq_index/PER_WORD]|= 1<<(seq_index%PER_WORD);
				j= block->start[seq_index+1];
			}
			else j++;     
		}
	}

	mask= 0;

	for (k=0;k<=(max_flex-flex);k+=max(1,(max_flex-flex)))
	{
		for (l=0;l<=k;l++)
		{
			act_d= min_dist-l;

			if ((act_d<min_d) || (act_d+flex+k>max_d))
				continue;

			for (num_seq=m=0;m<Length_Seq_Hit_Vector;m++) 
			{
            seq_mask= 0;
            for (n=act_d;n<=(act_d+k+flex);n++) 
				{
					assert (n-min_d>=0 && n-min_d <=2*max_flex);
               seq_mask|= Seq_Hit_Vector[n-min_d][m];
				}
            num_seq+= NUM_ONES32(seq_mask);
         }
			if (num_seq>=B_Options->Min_Nr_Seqs_Matching)
				mask|= (1<<(k+flex+16) | 1<<l);
		}
	}
		
	return mask;
}
	

/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
