/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef HITDEF
#define HITDEF

typedef struct t_hit_entry {
   void               *Pat;
   float               Information;
   struct t_hit_entry *Next;
   struct t_hit_entry *Prev;
} t_Hit_Entry;

typedef struct
{
   int          Num;
   int          Max_Num;
   t_Hit_Entry *First;
   t_Hit_Entry *Last;
	float        Min_Significance;
} t_Hits;

#ifdef HIT_MAIN
#define HIT_DEF
#else
#define HIT_DEF extern
#endif

HIT_DEF t_Hits *Init_Hit_List(int,float);
HIT_DEF int Hit_List_Will_Insert_Entry(t_Hits *, float );
HIT_DEF void Hit_List_Insert_Entry(t_Hits*,void *,float,int);
HIT_DEF int Pattern_Already_In_List(t_Hits *, void *);
HIT_DEF void Hit_List_Free(t_Hits *);


#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
