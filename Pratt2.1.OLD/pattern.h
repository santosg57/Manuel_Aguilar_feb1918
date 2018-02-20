/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef PATDEF
#define PATDEF

#define NFA_MAX_FLEX 20
#define MAX_WORD_NFA 20
#define PER_WORD_NFA 16

/** The limits correspond to a maximum pattern length of 320 positions */

typedef struct {
   int max_shortcut;
   int length;
   unsigned int char_vector[26];
   unsigned int shift_vector[NFA_MAX_FLEX];
} NFA_Pat1;

typedef struct {
   int max_shortcut;
   int length;
	int nr_words;
   unsigned int char_vector[26][MAX_WORD_NFA];
   unsigned int shift_vector[NFA_MAX_FLEX][MAX_WORD_NFA];
} NFA_Pat_longer;

typedef enum {shortnfa, longnfa} t_NFA_type;

typedef struct {
	t_NFA_type   type;
	union {
		NFA_Pat1        *pat1;
		NFA_Pat_longer  *patlong;
	} NFA;
} t_NFA;


#define PDEBUG 4567
typedef struct t_pattern
{
   int                debug;
   int                length;
	int                num_flex;
   int                num_comp;
   int               *positions;
   unsigned int      *spec_pos;
   int               *gaps;
   int               *flexi;
   int                num_sub;
   int                num_sub_new;
   int               *lengths;
   unsigned int     **bit_vecs;
   float              info;
	float              info_MDL;
	int                last;
	int                cardinality;
	int                num_seq;
} t_Pattern;




#define MAX_MATCHING_STRINGS 100

typedef struct t_pat_info {
   int                 num_comp;
   int                *positions;
   unsigned int       *spec_pos;
   int                *gaps;
   int                *flexi;
   int                *nflexi;
   int                 length;
	t_NFA              *NFA;
	int                 nr_match;
	int                 num_seq;
	int                 nr_occ;
	int                 nr_match_swiss;
	int                 nr_occ_swiss;
   unsigned int       *bit_vector;
   int                 num_sub;
   unsigned int      **bit_vecs;
	char              **matching_strings;
	int                *seq_indexes;
	int                 nr_matching_strings;
	int                 max_nr_matching_strings;
	struct t_pat_info  *refined_from;
	float               specificity;
	float               sensitivity;
	float               PPV;
	float               correlation;
	float               information;
	float               info_MDL;
	float               seq_hit_divergence;
	float               seq_hit_div_mst;
	float               fitness;
	unsigned int       *used;
	int                 special_flag;
} t_Pat_Info;

typedef struct {
	float               score;
	int                 num_comp;
   int                *positions;
   int                *gaps;
   int                *flexi;
   int                 length;
	int                 num_sub;
	int                 num_flex;
} t_Pat_Lead_Info;


#ifdef PAT_MAIN
#define PAT_DEF 
#else
#define PAT_DEF extern
#endif




PAT_DEF void Pat_Print_gcg(FILE*,t_Pat_Info *);
PAT_DEF void Pat_Print(FILE* ,t_Pat_Info *);
PAT_DEF void Pat_Print_Site(FILE* ,t_Pat_Info *);
PAT_DEF int Pat_Copy (t_Pattern*, t_Pattern*);
PAT_DEF void Pat_Check(t_Pattern*);
PAT_DEF t_Pat_Info* Pat_Info_From_Pattern(t_Pattern*,int);
PAT_DEF t_Pattern** Patterns_From_Alignment(t_alignment*,int*);

PAT_DEF void Pat_Info_Free(t_Pat_Info*);
PAT_DEF int Pattern_Equal_Pat_Info(t_Pattern*,t_Pat_Info*);
PAT_DEF int Pat_Info_Match_Same_Segments(t_Pat_Info*,t_Pat_Info*);

PAT_DEF int Pat_Info_From_String(char *, t_Pat_Info *);

PAT_DEF t_Pat_Lead_Info* Lead_Info_From_Pat(t_Pattern*,t_Pat_Lead_Info*);


PAT_DEF int Pattern_Num_X(t_Pattern*);
PAT_DEF int Pattern_Num_Char(t_Pattern*);
PAT_DEF int Pat_Info_Num_X(t_Pat_Info*);
PAT_DEF int Pat_Info_Num_Char(t_Pat_Info*);

PAT_DEF void Special_Print_Pattern_Matches(FILE* ,void* ,t_sequence** , int , t_Pat_Info* );


/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

#define G_MAX_COMP 101
#define G_MAX_FLEX 10

struct t_pat_edge;

typedef struct t_pat_node {
	int                 length_longest_path;
	float               max_score;
	float               Max_Score_Obtained;
   int                 Analysed;
	t_Pat_Lead_Info    *Best_Pat;
	int                 symbol;
	int                *pos_alignment;
	struct t_pat_edge **edges;
	int                 nr_edges;
	int                 max_nr_edges;
	float               Max_Score[G_MAX_COMP][G_MAX_FLEX];
	int                 seq_index;
} t_Pat_Node;
/*
	Node->Max_Score[c][n] is the maximum pattern score when
	starting from Node and allowing max c+1 components and max n
	flexibilities
*/

typedef struct t_pat_edge {
	unsigned int        possible[10][10];
	int                 min_aa;
	int                 max_aa;
	struct t_pat_node  *Node;
} t_Pat_Edge;

typedef struct {
	t_alignment        *alignment;
	int                 nr_seq_alignment;
	int                 max_nr_nodes;
	int                 nr_nodes;
	t_Pat_Node        **nodes;
} t_Pat_Graph;

struct t_block_struct;

PAT_DEF t_Pat_Graph* Pat_Graph_From_Alignment(t_alignment*,struct t_block_struct*);


#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
