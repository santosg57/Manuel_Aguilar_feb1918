/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#ifndef MENUDEF
#define MENUDEF

typedef enum {off,complete,start} t_restrict;

typedef struct
{
	int          Nr_Symbol_Block;
	int          Nr_Symbol_Search1;
	int          Full_Refinement;
	float        Min_Information;
	int          Gap_Restr;
	int          Max_Num_Flex;
	int          Max_Flex;
	int          Max_Flex_Prod;
	int          Max_Gap;
	int          Max_Length_menu;
	int          Max_Length;
	int          Max_Num_Comp;
	int          Min_Nr_Seqs_Matching;
	char        *Filename;
	int          Length_Hit_List;
	int          Sloppy;
	int          Number_Of_Occ_Output;
	int          Refinement;
	int          Refine_Generalise;
	int          Use_Short_Sequence_To_Guide_Search;
	int          Alignment_flag;
	char        *Al_Filename;
	int          Diagnostic;
	char        *Swiss_Flat_File;
   int          Quotient;
	int          Quotient_Inc;
	int          Print_Seqs_With_Motifs;
	int          Print_Ratio;
	int          Print_Vertical;
	int          Tree_Input;
	char        *Tree_Filename;
	int          Dist_Input;
	char        *Dist_Filename;
	int          Query_Mode;
	char        *Query_File_Name;
	int          MDL_flag;
	float        MDL_constant0;
	float        MDL_constant1;
	float        MDL_constant2;
	float        MDL_constant3;
	char        *filename;
	int          Show_Menu;
	int          Prosite_Style;
	t_restrict   Restrictions_Input;
	char        *Restrictions_File;
	int          Automatic;
	int          WWW;
	int          Max_Time;
	int          Input_Symbol_File;
	char         Symbol_File[100];
} t_Block_Options;


#ifdef MENU_MAIN
#define MENU_DEF
#else
#define MENU_DEF extern
#endif

MENU_DEF t_Block_Options* Menu_Set_Options(int,char*,t_Block_Options*);
MENU_DEF t_Block_Options* Options_From_File(char*);
MENU_DEF t_Block_Options* Options_From_Argumentline(int,char*,int,char**);
MENU_DEF void Options_Write_To_File(FILE*,t_Block_Options *);
MENU_DEF int Menu_Summarise_Options(t_Block_Options*);
MENU_DEF int OK_Parameter_Setting(t_Block_Options*);

#endif
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
