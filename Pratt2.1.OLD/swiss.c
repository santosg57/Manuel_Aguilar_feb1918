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
#include "block.h"
#include "menu.h"
#include "scan.h"
#define SWISS_MAIN
#include "swiss.h"

#define MAX_L 10000

extern t_Block_Options *B_Options;

int Find_Hits_In_SWISS_PROT(t_Hits* hits, char *filename)
/******************************************************************/
/* Reads and scans SWISS-PROT in normal format                    */
/******************************************************************/
{
	int i;
   t_Hit_Entry *Entry;
	char buf[100];
	char name[50];
	int line;
	int length;
   char seq_string[MAX_L];
   FILE *file;
   char *ok;
   char *seq;
   t_Pat_Info *pat;
	int nr_entries;

	printf ("Finding Hits in SWISS-PROT\n");

	nr_entries= 0;
   for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next) {
      if (!Entry->Pat)
         continue;
      pat= Entry->Pat;
      if ((pat->NFA= NFA_From_Pattern(Entry->Pat))==NULL)
         printf ("rotten pattern, could not make NFA!!\n");
      pat->nr_occ_swiss= pat->nr_match_swiss= 0;
		nr_entries++;
   }

	if (nr_entries==0)
		return 1;

	nr_entries= 0;

	printf ("scanning %s\n",filename);
	printf (" -- one dot per 1000 sequences scanned\n");
   file= fopen(filename,"r");

   if (!file) {
      fprintf (stderr,"Hit_Find_Hits_In_SWISS_PROT: cannot open file %s\n",filename);
		exit (1);
   }

	line= 0;
   ok=filename;
   while (ok!=NULL) {
      while ((ok=fgets(buf,100,file)) && (buf[0]!='I' || buf[1]!='D'))
         line++;
      if (!ok)
         continue;
      sscanf(buf+5,"%s",name);
      while ( (ok=fgets(buf,100,file)) && (buf[0]!='S' || buf[1]!='Q'))
         line++;
      if (!ok)
         continue;
		length= 0;
      while ( (ok=fgets(buf,100,file)) && (buf[0]!='/' || buf[1]!='/'))
      {
         line++;
         for (i=0;i<(int)strlen(buf);i++)
            if (isupper(buf[i]))
               seq_string[length++]= buf[i];
      }
      if (!ok) {
         printf ("sequence %s not terminated\n",name);
			printf ("line %d in %s\n",line,filename);
         return 0;
      }

		seq_string[length]='\0';
		nr_entries++;

      for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next) {
         pat= Entry->Pat;
         if (pat)
         {
            assert (pat->NFA);
            Num_Match_String_NFA(seq_string,pat,&(pat->nr_occ_swiss),&(pat->nr_match_swiss),B_Options->Refinement);
         }
      }
		if ((nr_entries+1)%1000==0)
			printf (".\n");
   }
	printf ("File %s with %d entries scanned\n",filename,nr_entries);
	if (nr_entries==0) {
		fprintf (stderr,"File %s had no SWISS-PROT format entries\n");
		fprintf (stderr,"No discriminatory analysis can be done\n");
		exit(1);
	}
   fclose (file);
   return nr_entries;
}



int Find_Hits_In_SWISS_PROT_Special_Format(t_Hits* hits, char *filename)
/******************************************************************/
/* Assumes one entry on each line with the entry name first, then 
   two spaces and then the amino acid sequence                    */
/******************************************************************/
{
   t_Hit_Entry *Entry;
   char seq_string[MAX_L];
   FILE *file;
   char *ok;
   char *seq;
   t_Pat_Info *pat;
	int nr_entries;

	printf ("Finding Hits in SWISS-PROT\n");

	nr_entries= 0;
   for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next) {
      if (!Entry->Pat)
         continue;
      pat= Entry->Pat;
      if ((pat->NFA= NFA_From_Pattern(Entry->Pat))==NULL)
         printf ("rotten pattern, could not make NFA!!\n");
      pat->nr_occ_swiss= pat->nr_match_swiss= 0;
		nr_entries++;
   }

	if (nr_entries==0)
		return 1;


	nr_entries= 0;

	printf ("scanning\n");
   file= fopen(filename,"r");

   if (!file) {
      fprintf (stderr,"Hit_Find_Hits_In_SWISS_PROT: cannot open file %s\n",filename);
      return 0;
   }

   ok=fgets(seq_string,MAX_L,file);
   while (ok!=NULL) {
      seq= strstr(seq_string,"  ");
      if (seq==NULL)
      {
         printf ("Empty sequence:%s\n",seq_string);
         ok=fgets(seq_string,MAX_L,file);
         continue;
      }
		nr_entries++;
      for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next) {
         pat= Entry->Pat;
         if (pat)
         {
            assert (pat->NFA);
            Num_Match_String_NFA(seq+2,pat,&(pat->nr_occ_swiss),&(pat->nr_match_swiss),1);
         }
      }
      ok=fgets(seq_string,MAX_L,file);
   }
   fclose (file);
   return nr_entries;
}

int Scan_Swiss_Prot_Refined_Pattern(t_Pat_Info* pat)
/******************************************************************/
{
	t_Pat_Info *Original_Pat;
	int i,j,k;
	int nr1,nr2;
	int last;

	if (!pat)
	{
		printf ("Scan_Swiss_Prot_Refined_Pattern: No pattern!!\n");
		return 0;
	}

   if ((pat->NFA= NFA_From_Pattern(pat))==NULL)
      printf ("rotten pattern, could not make NFA!!\n");

   pat->nr_occ_swiss= pat->nr_match_swiss= 0;

	Original_Pat= pat->refined_from;
	last= (-1);

	for (i=0;i<Original_Pat->nr_matching_strings;i++)
	{
		nr1=nr2= 0;
		Num_Match_String_NFA(Original_Pat->matching_strings[i],pat,&nr1,&nr2,0);
		if (nr2>0 && Original_Pat->seq_indexes[i]!=last) {
			pat->nr_match_swiss++;
			last= Original_Pat->seq_indexes[i];
		}
		pat->nr_occ_swiss+= nr2;
	}
	return 1;
}



static int Find_Hits_Refined_Patterns(t_Hits* hits)
/******************************************************************/
{
	t_Hit_Entry *Entry;

   for (Entry= hits->First; (Entry!=NULL); Entry= Entry->Next) {
		Scan_Swiss_Prot_Refined_Pattern(Entry->Pat);
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
