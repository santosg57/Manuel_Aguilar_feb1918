/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef BLOCKDEF
#define BLOCKDEF

#define PER_WORD 32


#define MAX_NR_SETS 1000

typedef struct t_block_struct
{
   int            len;
   int            alfa;
   int           *seq;
   int           *start;
   int           *length;
   unsigned int **vectors;
} t_block;

#ifdef BLOCK_MAIN
#define BLOCK_DEF
#else
#define BLOCK_DEF extern
#endif

#define MAX_NR_GEN 10
typedef struct 
{
	int           nr_symbols;
	int           symbols[MAX_NR_GEN];
} t_Generalize_Symbol_List;

BLOCK_DEF int Block_Init_Sets(char*);
BLOCK_DEF t_block* Block_Init(int, t_sequence *Seqs[]);
BLOCK_DEF void Block_Print_Hits(FILE*,t_block*,t_Pat_Info*);
BLOCK_DEF void Block_Print_Hits_To_Sequences(char,t_block*,t_Pat_Info*);
BLOCK_DEF void Statistics_Refined_Pattern(t_block*,t_Pat_Info*);
BLOCK_DEF void Hit_List_Find_Statistics(t_block*,t_Hits*);
BLOCK_DEF void Hit_Print_Hits(t_block*,t_Hits*,char*,char*,int );
BLOCK_DEF void Matching_Prepare(int,int);
BLOCK_DEF unsigned int Nr_Match(int,int,int,int,int,t_block *);

BLOCK_DEF int Pat_Info_Divergence(t_Pat_Info *,t_block *);

BLOCK_DEF void Hit_Print_mdl_Hits(t_block*,t_Hits*,FILE*);
#endif

#define LOOK_MASK 65535
#define LOOK_UP_ONES 65536
#define NUM_ONES32(n) (Num_Ones[(n)&LOOK_MASK]+Num_Ones[(n)>>16])





/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
