/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef SCANDEF
#define SCANDEF

#define WORD 32
#define MAX_WORD 8

#define MAX_FLEX 10
#define MAX_LEN_PAT 256

#define MAX_SYMBOL 20

#ifdef SCAN_MAIN
#define SCAN_DEF
#else
#define SCAN_DEF extern
#endif

SCAN_DEF t_NFA* NFA_From_Pattern( t_Pat_Info*);
SCAN_DEF int Num_Match_String_NFA(char*,t_Pat_Info*,int*,int*,int);

#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
