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
#include "al.h"
#include "pattern.h"
#define HIT_MAIN
#include "hit.h"

extern int Nr_Patterns_Evaluated;
extern void flush_status(t_Hits* );


extern int nr_now, nr_total;

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
#define fidentical(a,b) (fabs((a)-(b))<0.01)
#define fbigger(a,b) (((a)-(b))>0.01)

/*
#define DEBUG
*/

#ifdef DEBUG
#define CHECK(a) assert(a)
#else
#define CHECK(a)
#endif


float Hit_Needed;


t_Hits *Init_Hit_List(int max_num,float min_sig)
/******************************************************************/
{
	t_Hits *Hits;
	t_Hit_Entry *Entry;

	Hits= (t_Hits*)malloc(sizeof(t_Hits));
	assert (Hits);
	Hits->First= (t_Hit_Entry*)malloc(sizeof(t_Hit_Entry));
	assert (Hits->First);
	Hits->First->Pat= NULL;
	Hits->First->Information= 0.0;
	Hits->First->Next= NULL;
	Hits->First->Prev= NULL;
	Hits->Last= Hits->First;
	Hits->Max_Num= max_num;
	Hits->Num= 1;
	Hits->Min_Significance= min_sig;

	Hit_Needed= min_sig;

	return Hits;
}

int Pattern_Already_In_List(t_Hits *Hits, void *Pattern)
/******************************************************************/
{
	t_Hit_Entry *temp;

   temp= Hits->First;
   while (temp && fbigger(temp->Information,((t_Pattern*)Pattern)->info)) {
      temp= temp->Next;
   }
	while (temp && fidentical(((t_Pattern*)Pattern)->info,temp->Information)) {
		if (Pattern_Equal_Pat_Info(Pattern,temp->Pat))
			return 1;
		temp= temp->Next;
	}
	return 0;
}

int Hit_List_Will_Insert_Entry(t_Hits *Hits, float significance)
/******************************************************************/
{
	if ( (significance < Hits->Min_Significance) ||
		  ((significance <= Hits->Last->Information)&&(Hits->Num==Hits->Max_Num)))
		return 0;
	else
		return 1;
}

void Hit_List_Insert_Entry(t_Hits *Hits,void *Pattern, float significance,int refined_flag)
/******************************************************************/
{
	t_Hit_Entry *New_Entry;
	t_Hit_Entry *temp;

	/*
	Pat_Print(stdout,Pattern); 
	printf ("\n");
	*/

	if ( (significance < Hits->Min_Significance) ||
		  ((significance <= Hits->Last->Information)&&(Hits->Num==Hits->Max_Num)))
		return;

	New_Entry= (t_Hit_Entry*)malloc(sizeof(t_Hit_Entry));
	New_Entry->Information= significance;
	New_Entry->Pat= Pattern;

   if (New_Entry->Information>Hits->First->Information) {
		/* The new entry goes on top of the list */
      Hits->First->Prev= New_Entry;
      New_Entry->Next= Hits->First;
      New_Entry->Prev= NULL;
      Hits->First= New_Entry;
		/* Remove less significant patterns matching the same segments */
		temp=New_Entry->Next;
		while (temp && temp->Pat) {
			if (Pat_Info_Match_Same_Segments(New_Entry->Pat,temp->Pat)) {
				/* Remove temp from list */
				temp->Prev->Next= temp->Next;
				if (temp->Next==NULL) {
					Hits->Last= temp->Prev;
				}
				else {
					temp->Next->Prev= temp->Prev;
				}
				Hits->Num--;
				/*
				Pat_Info_Free(temp->Pat);
				*/
			}
			temp= temp->Next;
		}	
   }
   else
   {
		/* Find correct location of new entry, do not insert it if there is
			a more significant pattern matching the same segments */
      temp= Hits->First;
      while (temp->Next && New_Entry->Information<=temp->Next->Information) {
			if (Pat_Info_Match_Same_Segments(temp->Pat,New_Entry->Pat))
				return;
         temp= temp->Next;
      }
		/* Do not include new entry in list if a more significant pattern
			match the same segments */
		if (Pat_Info_Match_Same_Segments(temp->Pat,New_Entry->Pat))
			return;

		if (temp->Next) {
      	temp->Next->Prev= New_Entry;
      	New_Entry->Next= temp->Next;
      	New_Entry->Prev= temp;
      	temp->Next= New_Entry;
		}
		else {
		/* should this ever happen??? */
			temp->Next= New_Entry;
			New_Entry->Prev=temp;
			New_Entry->Next= NULL;
			Hits->Last= New_Entry;
		}
		/* Remove less significant patterns matching the same segments */
		temp=New_Entry->Next;
		while (temp && temp->Pat) {
			if (Pat_Info_Match_Same_Segments(New_Entry->Pat,temp->Pat)) {
				/* Remove temp from list */
				temp->Prev->Next= temp->Next;
				if (temp->Next==NULL) {
					Hits->Last= temp->Prev;
				}
				else {
					temp->Next->Prev= temp->Prev;
				}
				/*
				Pat_Info_Free(temp->Pat);
				*/
				Hits->Num--;
			}
			temp= temp->Next;
		}	
   }
   if (Hits->Num < Hits->Max_Num)
      Hits->Num++;
   else {
      temp= Hits->Last;
		/*
		Pat_Info_Free(temp->Pat);
		*/
      Hits->Last= temp->Prev;
      Hits->Last->Next= NULL;
   }
	if (!refined_flag)
		flush_status(Hits);

	if (Hits->Num>=Hits->Max_Num)
		Hit_Needed= Hits->Last->Information;
}


void Hit_List_Free(t_Hits *Hits)
/******************************************************************/
{
	t_Hit_Entry *temp,*temp2;

	if (Hits==NULL)
		return;

   temp= Hits->First;
   while (temp) {
		if (temp->Pat)
			Pat_Info_Free(temp->Pat);
		temp2= temp;
      temp= temp->Next;
		free (temp2);
   }
	free (Hits);
}
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
