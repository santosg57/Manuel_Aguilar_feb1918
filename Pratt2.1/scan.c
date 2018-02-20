/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "hit.h"
#include "sequence.h"
#include "al.h"
#include "pattern.h"
#define SCAN_MAIN
#include "scan.h"


#define max(a,b) ((a)>(b)?(a):(b))

extern unsigned int Member_In_Set[];

#define MAX_NAME 50
#define MAX_L 10000         /* maximum length one sequence */

/*********** LOCAL ROUTINES ***********************************************/
static int Num_Match_String_1(char*,t_Pat_Info*,int*,int*,int);
static int Num_Match_String_Long(char*,t_Pat_Info*,int*,int*,int);
static void Add_String_To_Pat_Info(t_Pat_Info*,char*, char*,int);


int Num_Match_String_NFA(char *string, t_Pat_Info* pat, int *nr_match, int *nr_seq_matched, int flag)
/**************************************************************************/
{
	if (pat->NFA->type==shortnfa)
		return Num_Match_String_1(string, pat,nr_match,nr_seq_matched,flag);
	else
		return Num_Match_String_Long(string,pat,nr_match,nr_seq_matched,flag);
}

t_NFA* NFA_From_Pattern( t_Pat_Info* pattern)
/**************************************************************************/
{
	int h,i,j,k;
	unsigned int NFA_pattern[MAX_LEN_PAT];
	int NFA_pattern_shortcut[MAX_LEN_PAT];
	int nr_pos;
	int nr_words;
	int symbol;
	unsigned int mask;
	t_NFA* temp;
	int upper;
	int max_short;

	max_short= 0;
	temp= malloc(sizeof(t_NFA));

	for (nr_pos=i=0;i<pattern->num_comp;i++)
	{
		NFA_pattern[nr_pos]= 0;
		if (pattern->positions[i]>=0) {
			symbol= pattern->positions[i];
			for (j=0;j<26;j++) {
	         if (Member_In_Set[symbol]&(1<<j)) {
					NFA_pattern[nr_pos]|= (1<<j);
            }
			}
		}
		else
		{
			mask=pattern->spec_pos[i];
			for (j=0;j<20;j++) {
				if (mask&(1<<j)) {
					for (k=0;k<26;k++) {
						if (Member_In_Set[j]&(1<<k)) {
							NFA_pattern[nr_pos]|=(1<<k);
						}
					}
				}
			}
		}
		if (i==pattern->num_comp-1)
			NFA_pattern_shortcut[nr_pos++]= 0;
		else 
		{
			if (pattern->gaps[i+1]==0)
				NFA_pattern_shortcut[nr_pos++]= pattern->flexi[i+1];
			else 
				NFA_pattern_shortcut[nr_pos++]= 0;

			upper= pattern->gaps[i+1]+pattern->flexi[i+1];
			for (j= 1;j<=upper;j++)
			{
				NFA_pattern[nr_pos]= ~0;
				if (j>=pattern->gaps[i+1])
					NFA_pattern_shortcut[nr_pos]= upper-j;
				else			
					NFA_pattern_shortcut[nr_pos]= 0;
				nr_pos++;
			}
		}
	}

	for (i=0;i<nr_pos;i++) {
		max_short= (max_short>=NFA_pattern_shortcut[i] ? max_short : NFA_pattern_shortcut[i]);
	}
	assert (max_short<NFA_MAX_FLEX);

	nr_words= nr_pos/WORD + (nr_pos%WORD?1:0);

	if (nr_words==1) {
		NFA_Pat1* NFA;
		temp->type= shortnfa;
		temp->NFA.pat1= NFA=  malloc(sizeof(NFA_Pat1));
		NFA->max_shortcut= max_short;

		NFA->length= nr_pos;
		for (i=0;i<26;i++){
			NFA->char_vector[i]= 0;
			for (j=0;j<nr_pos;j++)
				if (NFA_pattern[j]&(1<<i))
					NFA->char_vector[i]|= 1<<j;
		}
		for (i=0;i<=NFA->max_shortcut;i++) {
			NFA->shift_vector[i]= 0;
			for (j=0;j<nr_pos;j++)
				if (NFA_pattern_shortcut[j]>=i)
					NFA->shift_vector[i]|= 1<<j;
		}
	}
	else {
		NFA_Pat_longer* NFAL;
		temp->type= longnfa;
		temp->NFA.patlong= NFAL= malloc(sizeof(NFA_Pat_longer));
		NFAL->max_shortcut= max_short;

		NFAL->nr_words= nr_words= nr_pos/PER_WORD_NFA + (nr_pos%PER_WORD_NFA?1:0);
		assert (nr_words<MAX_WORD_NFA);
		NFAL->length= nr_pos;
		for (i=0;i<26;i++){
			for (j=0;j<nr_words;j++)
				NFAL->char_vector[i][j]= 0;
			for (j=0;j<nr_pos;j++)
				if (NFA_pattern[j]&(1<<i))
					NFAL->char_vector[i][j/PER_WORD_NFA]|= 1<<(j%PER_WORD_NFA);
		}
		for (i=0;i<=NFAL->max_shortcut;i++) {
			for (j=0;j<nr_words;j++)
				NFAL->shift_vector[i][j]= 0;
			for (j=0;j<nr_pos;j++)
				if (NFA_pattern_shortcut[j]>=i)
					NFAL->shift_vector[i][j/PER_WORD_NFA]|= 1<<(j%PER_WORD_NFA);
		}
	}

	return temp;
}


		


static int Num_Match_String_1(char *string, t_Pat_Info *pat, int *nr_match, int *nr_seq_matched,int flag)
/**************************************************************************/
{
	NFA_Pat1* NFA;
	unsigned int match[2];
	unsigned int temp;
	int h,i,j;
	char *cpt;
	int si;
	int current;
	unsigned int *cm;
	unsigned int accepting;
	int this_seq_matched_already;
	char match_string[100];
	static int index= 0;

	index++;

	NFA= pat->NFA->NFA.pat1;

	this_seq_matched_already= 0;

	accepting= 1<<(NFA->length-1);

	for (i=0;i<2;i++) 
		match[i]= 0;
	current= 0;
	for (cpt= string;(*cpt);cpt++) {
		cm=&(match[current]);
		si= (*cpt)-'A';
		temp= match[1-current];
		*cm= temp<<1 | 1;
		for (i=1;i<=NFA->max_shortcut;i++)
			*cm|= (temp & NFA->shift_vector[i])<<(i+1);
		*cm&= NFA->char_vector[si];
		if (*cm&accepting) {
			if (flag)
				Add_String_To_Pat_Info(pat, cpt, string,index);
				
		/*
			strncpy(match_string,cpt-30,31);
			printf ("%s match\n",match_string);
		*/
			(*nr_match)++;
			if (!this_seq_matched_already)
				(*nr_seq_matched)++;
			this_seq_matched_already= 1;
		}
		current= 1-current;	
	}
}

static int Num_Match_String_Long(char *string, t_Pat_Info* pat, int *nr_match, int *nr_seq_matched,int flag)
/**************************************************************************/
{
	NFA_Pat_longer* NFA;
	unsigned int match[2][MAX_WORD_NFA];
	unsigned int temp;
	int h,i,j,k;
	char *cpt;
	int si;
	int current;
	unsigned int *cm;
	unsigned int accepting;
	int this_seq_matched_already;
	char match_string[100];
	unsigned int overflow;
	static int index= 0;

	index++;

	NFA= pat->NFA->NFA.patlong;

	this_seq_matched_already= 0;

	accepting= 1<<((NFA->length-1)%PER_WORD_NFA);

	for (i=0;i<2;i++) 
		for (j=0;j<NFA->nr_words;j++)
			match[i][j]= 0;
	current= 0;
	for (cpt= string;(*cpt);cpt++) {
		cm=&(match[current][0]);
		si= (*cpt)-'A';
		overflow= 1;
		for (i=0;i<NFA->nr_words;i++)
		{
			temp= match[1-current][i];
			*cm= temp<<1 | overflow;
			for (j=1;j<=NFA->max_shortcut;j++)
				*cm|= (temp & NFA->shift_vector[j][i])<<(j+1);
			overflow= *cm >> PER_WORD_NFA;
			*cm&= NFA->char_vector[si][i];
			if (i==NFA->nr_words-1 && *cm&accepting) {
				if (flag)
					Add_String_To_Pat_Info(pat, cpt, string,index);
			/*
				strncpy(match_string,cpt-NFA->length,NFA->length+1);
				printf ("%s match\n",match_string);
			*/
				(*nr_match)++;
				if (!this_seq_matched_already)
					(*nr_seq_matched)++;
				this_seq_matched_already= 1;
			}
			cm++;
		}
		current= 1-current;	
	} 
} 

static void Add_String_To_Pat_Info(t_Pat_Info* pat, char *end_of_match, char *string, int index) 
/**************************************************************************/ 
{ 
	int i,j; 
	char *temp; 
	
	if (pat->nr_matching_strings==pat->max_nr_matching_strings) 
	{ 
		pat->max_nr_matching_strings+=100;
		pat->matching_strings=realloc(pat->matching_strings,sizeof(char*)*pat->max_nr_matching_strings);
		pat->seq_indexes=realloc(pat->seq_indexes,sizeof(int)*pat->max_nr_matching_strings);
		assert (pat->matching_strings);
	}
	pat->seq_indexes[pat->nr_matching_strings]= index;
	temp= pat->matching_strings[pat->nr_matching_strings++]= calloc(strlen(string)+5,sizeof(char));
	assert (temp);
	strcpy(temp,string);
}
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
