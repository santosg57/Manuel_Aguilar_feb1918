/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef AL_DEF
#define AL_DEF

typedef struct
{
	int nr_seq;
	int nr_pos;
	int *seq_index; /* for alignments from concatenated sequences
							 - to avoid patterns crossing sequence boundaries */
	char **positions;
	char **names;
} t_alignment;

t_alignment* Alignment_From_File (char* );
t_alignment* Alignment_From_Shortest_Sequence(int, t_sequence**,int);
t_alignment* Al_From_Sequences(t_sequence **,int);
void Alignment_Print (t_alignment*);
void Alignment_Print_Seqs_Fasta(t_alignment*,char*);


#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
