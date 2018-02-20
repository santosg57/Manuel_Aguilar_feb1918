/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef SWISSDEF
#define SWISSDEF

#ifdef SWISS_MAIN
#define SWISS_DEF 
#else
#define SWISS_DEF extern
#endif

SWISS_DEF int Find_Hits_In_SWISS_PROT(t_Hits* , char *);
SWISS_DEF int Scan_Swiss_Prot_Refined_Pattern(t_Pat_Info*);

#endif



/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
