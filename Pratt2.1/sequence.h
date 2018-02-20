/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef SEQREADDEF
#define SEQREADDEF

#define MAX_RESTRICT 10

typedef enum {dna,rna,protein,unknown}  t_seq_type;

typedef struct 
{
	int          nr;
	char        *base;
	int         *freq;
	float        ratio[128];
	float        exp;
} t_seq_freq;

typedef struct
{
	int             length;
	unsigned char  *seq;
	unsigned char  *pot;
	short          *tab;
} t_seq_graph;

typedef struct
{
	int          index;
	char        *filename;
	t_seq_type   type;
	int          length;
	int          max_length;
	char        *seq;
	char        *out;
	t_seq_freq  *freq;
	t_seq_graph *graph;
	int          pat_position;
	int          tree_index;
	/* flag which is set to 1 if matches to this sequence is restricted */
	int          restricted;
	/* index of first and last possible starting point matching a pattern */
	int          first[MAX_RESTRICT];
	int          last[MAX_RESTRICT];  
} t_sequence;


#define seq_type(s) (s)->type
/*
#define seq_length(s) (s)->length
*/
#define seq_element(s,i) ( (s)->seq[i] )
#define seq_element2(s,i) ( ((i>=0) && (i<(s)->length)) ? ((s)->graph->seq[i]):0) 
#define seq_exp(s) ((s)->freq)->exp
#define seq_filename(s) ((s)->filename)


#ifdef SEQ_MAIN
#define SEQ_DEF
char DNA_alphabet[]= "ACTGN\0";
char RNA_alphabet[]= "ACUG\0";
char PROTEIN_alphabet[]= "ARNDCQEGHILKMFPSTWYVBZX\0";
int Collins[26][26];
unsigned char Trans_Table[128];
#else
#define SEQ_DEF extern
extern char DNA_alphabet[];
extern char RNA_alphabet[];
extern char PROTEIN_alphabet[];
extern int Collins[26][26];
extern unsigned char Trans_Table[];
#endif


SEQ_DEF t_seq_graph* seq_graph_from_seq(t_sequence*);
SEQ_DEF int Sequence_Initialize_Translation();
SEQ_DEF t_sequence *seq_init_sequence(int,t_seq_type,char*);
SEQ_DEF void seq_add_res(t_sequence*,char);
SEQ_DEF t_sequence *seq_from_GCG_file(char*);
SEQ_DEF void seq_print_seq(t_sequence*,int);
SEQ_DEF void seq_print_part(t_sequence*,int,int);
SEQ_DEF void seq_print_part_file(FILE*,t_sequence*,int,int);
SEQ_DEF float seq_dist (char* , int , char *, int );
#define STORAGE_ERROR 1
#define FILE_OPEN_ERROR 2
#define ROTTEN_ERROR 3
SEQ_DEF void seq_error_message (int ,char* , int );
SEQ_DEF void Print_Sequence_To_File(FILE**,t_sequence*);

SEQ_DEF t_sequence** seqs_from_SwissProt_file(char *,int *,t_seq_type );
SEQ_DEF t_sequence** seqs_from_SRS_file(char *,int *,t_seq_type );
SEQ_DEF t_sequence** seqs_from_Fasta_file(char *,int *,t_seq_type );
SEQ_DEF t_sequence* one_seq_from_Fasta_file(char *,t_seq_type );

SEQ_DEF void Print_Seqs_To_Fasta_File(t_sequence**,int,char*);

SEQ_DEF void Seqs_Prepare_For_Output(t_sequence** , int );
SEQ_DEF int Seq_Add_Symbol_To_Output(t_sequence* , int , char );
SEQ_DEF int Seqs_Print_Out(FILE **,t_sequence** , int );

SEQ_DEF int Sequences_Read_Restrictions_File(t_sequence**,int,char*);
SEQ_DEF void Print_Restrictions(FILE*,t_sequence**,int);

SEQ_DEF void Print_Seq_Lengths(FILE*,t_sequence**,int);

#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
