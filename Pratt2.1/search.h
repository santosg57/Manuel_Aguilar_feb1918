/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef SEARCHDEF
#define SEARCHDEF

typedef struct {
   t_Pat_Info         *Pat;
   int                 pos;
   int                 gaps;
   int                 flexi;
} t_Search;

#ifdef SEARCH_MAIN
#define SEARCH_DEF
#else
#define SEARCH_DEF extern
#endif

SEARCH_DEF void Find_Motif_Block(t_Block_Options*);
SEARCH_DEF t_Hits* Refine_Hits(t_block *,t_Hits*);
SEARCH_DEF int Hit_Find_Hits_In_SWISS_PROT(t_Hits*,char*);
SEARCH_DEF t_Hits *Find_Block(t_block*);

#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
