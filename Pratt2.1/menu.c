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
#define MENU_MAIN
#include "menu.h"
#include "help.h"
#include "values.h"

extern int nrseqs;
extern int Total_Seq_Length;

extern usage();

static int name_set= 0;


#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

static t_Block_Options* Option_Init(int nr_seqs,char *name)
/******************************************************************/
{
	t_Block_Options *temp;
	
	temp= (t_Block_Options*)malloc(sizeof(t_Block_Options));
	assert (temp);
	temp->Min_Nr_Seqs_Matching= nr_seqs;
	temp->Nr_Symbol_Block= 20;
	temp->Nr_Symbol_Search1= 20;
	temp->Refinement= 1;
	temp->Full_Refinement= 0;
	temp->Refine_Generalise= 0;
	temp->Min_Information= -10000.0;
	temp->Sloppy= 0;
	temp->Gap_Restr= 1;
	temp->Max_Num_Flex= 2;
	temp->Max_Flex= 2;
	temp->Max_Flex_Prod= 10;
	temp->Max_Gap= 5;
	temp->Max_Length_menu= 50;
	temp->Max_Num_Comp= 50;
	temp->Length_Hit_List= 50;
	temp->Number_Of_Occ_Output= 50;
	temp->Print_Ratio= 10;
	temp->Input_Symbol_File= 0;
	sprintf(temp->Symbol_File,"Pratt.sets");

	temp->Show_Menu= 0;
	
	temp->filename=calloc(100,sizeof(char));
	strcpy(temp->filename,name);

	temp->Query_Mode= 0;
	temp->Query_File_Name=calloc(100,sizeof(char));
	sprintf(temp->Query_File_Name,"query");

	temp->MDL_flag= 0;
	temp->MDL_constant0= 10.0;
	temp->MDL_constant1= 10.0;
	temp->MDL_constant2= 3.0;
	temp->MDL_constant3= 10.0;

	temp->Print_Seqs_With_Motifs= 1;
	temp->Filename= (char*)malloc(100*sizeof(char));
	assert (temp->Filename);
	sprintf (temp->Filename,"%s.%d.pat\0",name,temp->Min_Nr_Seqs_Matching);

	temp->Alignment_flag= 0;
	temp->Al_Filename= (char*)malloc(40*sizeof(char));
	temp->Use_Short_Sequence_To_Guide_Search= 1;
	sprintf (temp->Al_Filename,"%s.aln\0",name);
	temp->Diagnostic= 0;
	temp->Quotient= 3;
	temp->Quotient_Inc= 10000;

	temp->Prosite_Style= 1;
	temp->Restrictions_Input= off;

	temp->Tree_Input= 0;
	temp->Tree_Filename= calloc(40,sizeof(char));
	sprintf (temp->Tree_Filename,"%s.dnd",name);
	temp->Dist_Input= 0;
	temp->Dist_Filename= calloc(40,sizeof(char));
	sprintf (temp->Dist_Filename,"%s.dist",name);
	temp->Swiss_Flat_File= calloc(100,sizeof(char));
	sprintf(temp->Swiss_Flat_File,"sprot32.dat");
	temp->Restrictions_File= calloc(100,sizeof(char));
	sprintf(temp->Restrictions_File,"%s.restr",name);

	temp->Automatic= 0;
	temp->Max_Time= MAXINT;

	return temp;
}

static int Get_Real_Value(float lower, float upper,float *value)
/******************************************************************/
{
	char buf[100];
	printf ("New value (real): ");
	scanf ("%f", value);
	if (*value<lower || *value>upper) {
		printf ("value must be between %d and %d\n",lower,upper);
		printf ("press return to continue\n");
		fgets (buf,100,stdin);
		fgets (buf,100,stdin);
		return 0;
	}
	return 1;
}


static int Get_Int_Value(int lower, int upper,int *value)
/******************************************************************/
{
	char buf[100];
	printf ("New value: ");
	scanf ("%d", value);
	if (*value<lower || *value>upper) {
		printf ("value must be between %d and %d\n",lower,upper);
		printf ("press return to continue\n");
		fgets (buf,100,stdin);
		fgets (buf,100,stdin);
		return 0;
	}
	return 1;
}



t_Block_Options* Menu_Set_Options(int nr_seqs,char *name, t_Block_Options* options)
/******************************************************************/
{
	t_Block_Options *temp;
	char buf[100];
	char buf2[100];
	char buf3[100];
	int i;
	int par;
	float fpar;

	if (options) {
		temp= options;
	}
	else 
		temp= Option_Init(nr_seqs,name);
	
	buf[0]='\0';
	
	while (buf[0]!='X' && buf[0]!='x')
	{
		Options_Write_To_File(stdout,temp);
		printf ("\n");
		printf ("X: eXecute program\n");
		printf ("Q: Quit\n");
		printf ("\n");
		printf ("help: for on-line help\n");
		printf ("\n");

		printf ("Command: ");
		scanf ("%s",buf);
		for (i=0;i<strlen(buf);i++)
			if (islower(buf[i]))
				buf[i]= toupper(buf[i]);

		switch (buf[0]) {

			case 'B' : 
				switch (buf[1]) {
					case 'N': 
						printf ("New value: ");
						scanf ("%d", &par);
						if (par<0)
							printf ("Negative values not accepted\n");
						else {
							temp->Nr_Symbol_Block= par; 
						   temp->Nr_Symbol_Search1= par; break;
						}
						break;
					case 'I':
						temp->Input_Symbol_File= 1-temp->Input_Symbol_File; 
						break;
					case 'F':
						if (temp->Input_Symbol_File) {
							printf ("Symbol filename: ");
							scanf ("%s",temp->Symbol_File);
							break;
						}
					default:
						fgets (buf3,100,stdin);
						printf ("Unknown command: %s\n",buf);
						printf ("Press return to continue.\n");
						fgets (buf3,100,stdin);
				}
				break;


			case 'C': 
				if (buf[1]=='M') {
					if (!Get_Int_Value(2,nrseqs,&par))
						break;
					temp->Min_Nr_Seqs_Matching= par; 
					if (!name_set)
						sprintf (temp->Filename,"%s.%d.pat\0",name,temp->Min_Nr_Seqs_Matching);
				}
				else if (buf[1]=='%') {
					if (!Get_Real_Value(1.0,100.0,&fpar))
						break;
					temp->Min_Nr_Seqs_Matching=(int)ceil((max(0.0,min(100.0,(double)fpar))*(double)nrseqs)/(double)100.0);  
					if (!name_set)
					 	sprintf (temp->Filename,"%s.%d.pat\0",name,temp->Min_Nr_Seqs_Matching);
				}
				else {
					fgets (buf3,100,stdin);
					printf ("Not a legal command (CM and C% are legal commands)\n");
					printf ("Press return to continue.\n");
					fgets (buf3,100,stdin);
				}
				break;



			case 'E': 
				if (Get_Int_Value(0,10000,&par))
					temp->Quotient= par; 
				break;


			case 'F':
				if (buf[1]!='N' && buf[1]!='L' && buf[1]!='P') {
					printf ("Not a legal command (FN, FL and FP are legal commands)\n");
					fgets (buf3,100,stdin);
					printf ("Press return to continue.\n");
					fgets (buf3,100,stdin);
					break;
				}
				if (Get_Int_Value(0,20,&par)) {
					switch (buf[1]) {
						case 'L': temp->Max_Flex= par; break; 
						case 'N': temp->Max_Num_Flex= par; break;
						case 'P': temp->Max_Flex_Prod= par; break;
						default: 
							fprintf(stderr,"error in %s at line %d\n",__FILE__,__LINE__),
							fprintf(stderr,"Please report to inge@ii.uib.n\n");
							exit(1);
					}
				}
				break;


			case 'G' : 
				if (buf[1]!='F' || (temp->Alignment_flag==0 && temp->Query_Mode==0)) {
					fgets (buf3,100,stdin);
					if (sscanf(buf3,"%s",buf2)==1) 
						;
					else {
						printf ("[seq,al,query]:");
						scanf ("%s", buf2);
						fgets (buf3,100,stdin);
					}
					for (i=0;i<strlen(buf2);i++)
						if (islower(buf2[i]))
							buf2[i]= toupper(buf2[i]);
					if (strcmp("SEQ",buf2)==0) {
						temp->Alignment_flag=temp->Query_Mode= 0;
					}
					else if (strcmp("AL",buf2)==0) {
						temp->Alignment_flag= 1;
						temp->Query_Mode= 0;
					}
					else if (strcmp("QUERY",buf2)==0) {
						temp->Alignment_flag= 0;
						temp->Query_Mode= 1;
					}
					else {
						printf ("Unknown option (only [seq,al,query] are acceptable)\n");
						printf ("Press return to return to menu\n");
						fgets (buf3,100,stdin);
					}
				}
				else if (buf[1]=='F') {
					if (temp->Alignment_flag) {
						printf ("Alignment file name: ");
						scanf ("%s",temp->Al_Filename);
					}
					else {
						printf ("Query file name: ");
						scanf ("%s",temp->Query_File_Name);
					}
				}
				break;


			case 'H':
				fgets (buf3,100,stdin);
				if (sscanf(buf3,"%s",buf2)==1)
					;
				else {
					printf ("Specify menu option or 'help' for general help:");
					scanf ("%s", buf2);
					fgets (buf3,100,stdin);
				}
				for (i=0;i<strlen(buf2);i++)
					if (islower(buf2[i]))
						buf2[i]= toupper(buf2[i]);

					specific_help(buf2);

				break;


			case 'M' : 
				if (!temp->Print_Seqs_With_Motifs || (buf[1]!='R' && buf[1]!='V')) {
					temp->Print_Seqs_With_Motifs= 1-temp->Print_Seqs_With_Motifs; 
				}
				else if (buf[1]=='R') {	
					if (Get_Int_Value(1,100,&par))
						temp->Print_Ratio=par;
				}
				else if (buf[1]=='V') {	
					temp->Print_Vertical=1-temp->Print_Vertical;
				}
				break;


			case 'O': 
				switch (buf[1]) {
					case 'F':
						name_set= 1;
						printf ("File name: ");	
						scanf ("%s", temp->Filename);
						break;
					case 'P': 
						temp->Prosite_Style= 1-temp->Prosite_Style; 
						break;
					case 'N': case 'A':
						printf ("New value: ");
						scanf ("%d", &par);
						if (par<0)
							printf ("Negative values not accepted\n");
						if (buf[1]=='N')
							temp->Length_Hit_List= par; 
						else
							temp->Number_Of_Occ_Output= par; 
						break;
					default:
						printf ("Not a legal command (OF, OP, ON and OA are legal commands)\n");
						fgets (buf3,100,stdin);
						printf ("Press return to continue.\n");
						fgets (buf3,100,stdin);
						break;
					}
				break;


			case 'P': 
				if (buf[1]!='L' && buf[1]!='N' && buf[1]!='X' && buf[1]!='P' && (temp->Restrictions_Input==off || buf[1]!='F')) {
					printf ("Not a legal command (PL, PN, PX, PP %s are legal commands)\n",((temp->Restrictions_Input!=off)?"and PF":""));
					fgets (buf3,100,stdin);
					printf ("Press return to continue.\n");
					fgets (buf3,100,stdin);
					break;
				}
				if (buf[1]=='P') {
					fgets (buf3,100,stdin);
					if (sscanf(buf3,"%s",buf2)==1) 
						;
					else {
						printf ("[off,complete,start]:");
						scanf ("%s", buf2);
						fgets (buf3,100,stdin);
					}
					for (i=0;i<strlen(buf2);i++)
						if (islower(buf2[i]))
							buf2[i]= toupper(buf2[i]);
					if (strcmp("OFF",buf2)==0)
						temp->Restrictions_Input= off;
					else if (strcmp("COMPLETE",buf2)==0)
						temp->Restrictions_Input= complete;
					else if (strcmp("START",buf2)==0)
						temp->Restrictions_Input= start;
					else {
						printf ("Unknown option (only [off,complete,start] are acceptable)\n");
						printf ("Press return to return to menu\n");
						fgets (buf3,100,stdin);
					}
				}
				else if (buf[1]=='F' && temp->Restrictions_Input!= off) {
					printf ("Restrictions filename\n");
					scanf("%s",temp->Restrictions_File);
				}
				else {
					if (!Get_Int_Value(0,100,&par))
						break;
					switch (buf[1]) {
						case 'X': temp->Max_Gap= par; break;
						case 'L': temp->Max_Length_menu= par; break;
						case 'N': temp->Max_Num_Comp= par; break;
						default: 
							fprintf(stderr,"error in %s at line %d\n",__FILE__,__LINE__),
							fprintf(stderr,"Please report to inge@ii.uib.n\n");
							exit(1);
							break;
					}
				}
				break;


			case 'Q': 
				printf ("Ciao\n");
				exit (1);
				break;


			case 'R': 
				if (buf[1]=='G')
					temp->Refine_Generalise= 1-temp->Refine_Generalise; 
				else
					temp->Refinement= 1-temp->Refinement; 
				break;


			case '7': 
				printf ("Restrictions filename\n");
				scanf("%s",temp->Restrictions_File);
				break;

			case 'S':
				if (buf[1]!='F' || (temp->Tree_Input==0 && temp->Dist_Input==0)) {
					fgets (buf3,100,stdin);
					if (sscanf(buf3,"%s",buf2)==1) 
						;
					else {
						printf ("[info,mdl,tree,dist,ppv]:");
						scanf ("%s", buf2);
						fgets (buf3,100,stdin);
					}
					for (i=0;i<strlen(buf2);i++)
						if (islower(buf2[i]))
							buf2[i]= toupper(buf2[i]);

					if (strcmp("TREE",buf2)==0) {
						temp->Tree_Input=1; temp->Diagnostic=temp->Dist_Input= temp->MDL_flag= 0;
					}
					else if (strcmp("DIST",buf2)==0) {
						temp->Dist_Input=1; temp->Diagnostic=temp->Tree_Input= temp->MDL_flag= 0;
					}
					else if (strcmp("MDL",buf2)==0) {
						temp->MDL_flag=1; temp->Diagnostic=temp->Tree_Input= temp->Dist_Input= 0;
					}
					else if (strcmp("PPV",buf2)==0) {
						temp->Diagnostic=1; temp->MDL_flag= temp->Tree_Input= temp->Dist_Input= 0;
					}
					else if (strcmp("INFO",buf2)==0) {
						temp->MDL_flag= temp->Tree_Input= temp->Dist_Input=temp->Diagnostic= 0;
					}
					else {
						printf ("Illegal parameter\n");
					}
				}
				else if (buf[1]=='F') {
					if (temp->Tree_Input) {
						printf ("Tree input filename: ");
						scanf ("%s",temp->Tree_Filename);
					}
					if (temp->Dist_Input) {
						printf ("Dist file name: ");
						scanf("%s",temp->Dist_Filename);
					}
					if (temp->Diagnostic) {
						printf ("SWISS-PROT flat filename: ");
						scanf("%s",temp->Swiss_Flat_File);
					}
				}
				break;


/*
			case 'U': 
				if (!Get_Int_Value(0,10000,&par))
					temp->Quotient_Inc= par;break;
				break; 
*/


			case 'X': case 'x': printf ("OK\n\n");
				break;


			case 'Z': 
				if (temp->MDL_flag) {
					if ((int)strlen(buf)==1)
						temp->MDL_flag= 0;
					else {
						printf ("New (real) value:");
						scanf("%f",&fpar);
						if (buf[1]=='0')
							temp->MDL_constant0= fpar;
						else if (buf[1]=='1')
							temp->MDL_constant1= fpar;
						else if (buf[1]=='2')
							temp->MDL_constant2= fpar;
						else if (buf[1]=='3')
							temp->MDL_constant3= fpar;
						else 
							fprintf (stderr,"no MDL parameter assignment made\n");
					}
				}
				else
					temp->MDL_flag= 1;
				break; 


			default:
				fgets (buf3,100,stdin);
				printf ("Unknown command: %s\n",buf);
				printf ("Press return to continue.\n");
				fgets (buf3,100,stdin);
		}
	}
	if (temp->Max_Num_Flex==0) {
		temp->Max_Flex= 0;
		temp->Max_Flex_Prod= 1;
	}

	temp->Max_Length=temp->Max_Length_menu+1;

	return temp;
}

t_Block_Options* Options_From_Argumentline(int nr_seqs,char *name, int argc, char **argv)
/******************************************************************/
{
	t_Block_Options *temp;
	int i,j;
	char option;
	int value;
	float fvalue;
	int on;
	char inupper[100];
	int isnext;
	char next[100];

	temp= Option_Init(nrseqs,name);

	i=3;
	while (i<argc) {
		if (argv[i][0]!='-')
			usage();

		if (strlen(argv[i])>=100)
			usage();

		strcpy(inupper,&(argv[i][1]));
		for (j=0;j<strlen(inupper);j++)
			if (islower(inupper[j]))
				inupper[j]= toupper(inupper[j]);

		if (i+1>=argc)
			isnext= 0;
		else {
			isnext= 1;
			strcpy(next,argv[i+1]);
			for (j=0;j<strlen(next);j++) {
				if (islower(next[j]))
					next[j]= toupper(next[j]);
			}
		}

		/* FIRST THE OPTIONS THAT CAN ONLY BE SET FROM COMMAND-LINE */
		if (strcmp(inupper,"MENU")==0) {
			temp->Show_Menu= 1;
			i++;
		}
		else if (strcmp(inupper,"AUTO")==0) {
			temp->Automatic= 1;		
			i++;
		}
		else if (strcmp(inupper,"WWW")==0) {
			temp->WWW= 1;		
			i++;
		}
		else if (strcmp(inupper,"MAXTIME")==0) {
         if (!isnext){
            fprintf (stderr,"After option %s, a positive integer is needed\n",argv[i]);
            exit(1);
         }
         temp->Max_Time= atoi(argv[i+1]);
         i+=2;
      }
		/* THEN THE OPTIONS THAT ALSO CAN BE SET FROM THE MENU - IN ALPHABETICAL ORDER */
		else if (strcmp(inupper,"BN")==0) {
         if (!isnext){
            fprintf (stderr,"After option %s, a positive integer is needed\n",argv[i]);
            exit(1);
         }
         temp->Nr_Symbol_Block=temp->Nr_Symbol_Search1= atoi(argv[i+1]);
			if (temp->Nr_Symbol_Block<1 || temp->Nr_Symbol_Block>100){
            fprintf (stderr,"%d is not allowed value for option %s\n",temp->Nr_Symbol_Block,argv[i]);
            exit(1);
         }
			i+=2;
		}
		else if (strcmp(inupper,"BI")==0) {
			if (!isnext) {
            fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"ON")==0)
				temp->Input_Symbol_File= 1;
			else if (strcmp(next,"OFF")==0)
				temp->Input_Symbol_File= 0;
			else {
            fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"BF")==0) {
			if (!isnext) {
            fprintf (stderr,"After option %s, a filename is expected\n",argv[i]);
				exit(1);
			}
			strcpy(temp->Symbol_File,argv[i+1]);
			i+=2;
		}
		else if (strcmp(inupper,"CM")==0) {
			if (!isnext) {
            fprintf (stderr,"After option %s, an integer is expected\n",argv[i]);
				exit(1);
			}
			temp->Min_Nr_Seqs_Matching= atoi(argv[i+1]);
			if (/* (nr_seqs>1) && */ ((temp->Min_Nr_Seqs_Matching<2) || (temp->Min_Nr_Seqs_Matching>nr_seqs))){
            fprintf (stderr,"%d is illegal min nr. of seqs. to match a patern (option %s)\n",temp->Min_Nr_Seqs_Matching,argv[i]);
				exit(1);
			}
			if (!name_set)
				sprintf (temp->Filename,"%s.%d.pat\0",name,temp->Min_Nr_Seqs_Matching);
			i+=2;
		}
		else if (strcmp(inupper,"C%")==0) {
			double ratio;
			if (!isnext) {
            fprintf (stderr,"After option %s, a real is expected\n",argv[i]);
				exit(1);
			}
			ratio= atof(argv[i+1]);
			if (ratio<=0.0 || ratio>100.0) {
				fprintf(stderr,"%s is illegal value for parameter %s\n",argv[i+1],argv[i]);
				exit(1);
			}
			temp->Min_Nr_Seqs_Matching=(int)ceil((max(0.0,min(100.0,ratio))*(double)nr_seqs)/(double)100.0);
			if (!name_set)
				sprintf (temp->Filename,"%s.%d.pat\0",name,temp->Min_Nr_Seqs_Matching);
			i+=2;
		}
		else if (strcmp(inupper,"E")==0) {
			if (!isnext) {
            fprintf (stderr,"After option %s, a non-negative integer is expected\n",argv[i]);
				exit(1);
			}
			temp->Quotient= atoi(argv[i+1]);
			if (temp->Quotient<0) {
            fprintf (stderr,"After option %s, a non-negative integer is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if ((strcmp(inupper,"FN")==0)||(strcmp(inupper,"FL")==0)|| (strcmp(inupper,"FP")==0)) {
			int par;
			if (!isnext) {
            fprintf (stderr,"After option %s, a non-negative integer is expected\n",argv[i]);
				exit(1);
			}
			par= atoi(argv[i+1]);
			if (par<0) {
            fprintf (stderr,"After option %s, a non-negative integer is expected\n",argv[i]);
				exit(1);
			}
			switch(inupper[1]) {
				case 'L': temp->Max_Flex= par; break;
            case 'N': temp->Max_Num_Flex= par; break;
            case 'P': temp->Max_Flex_Prod= par; break;
            default: exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"G")==0) {
			if (!isnext) {
            fprintf (stderr,"After option %s, one of SEQ, AL, or QUERY is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"SEQ")==0) {
				temp->Alignment_flag=temp->Query_Mode= 0;
			}
			else if (strcmp(next,"AL")==0) {
				temp->Alignment_flag= 1;
            temp->Query_Mode= 0;
			}
			else if (strcmp(next,"QUERY")==0) {
				temp->Alignment_flag= 0;
            temp->Query_Mode= 1;
			}
			else {
            fprintf (stderr,"After option %s, one of SEQ, AL, or QUERY is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"GF")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, a filename is expected\n",argv[i]);
				exit(1);
			}
			strcpy(temp->Al_Filename,argv[i+1]);
			strcpy(temp->Query_File_Name,argv[i+1]);
			i+=2;
		}
		else if (strcmp(inupper,"H")==0) {
			specific_help("H");
			exit(1);
		}
		else if (strcmp(inupper,"M")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, on or off is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"ON")==0)
				temp->Print_Seqs_With_Motifs= 1;
			else if (strcmp(next,"OFF")==0)
				temp->Print_Seqs_With_Motifs= 0;
			else {
				fprintf (stderr,"After option %s, on or off is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"MR")==0) {
			int par;
			if (!isnext) {
				fprintf (stderr,"After option %s, a positive integer value is expected\n",argv[i]);
				exit(1);
			}
			par=atoi(argv[i+1]);
			if (par<1) {
				fprintf (stderr,"After option %s, a positive integer value is expected\n",argv[i]);
				exit(1);
			}
			temp->Print_Ratio=par;
			i+=2;
		}
		else if (strcmp(inupper,"MV")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"ON")==0)
				temp->Print_Vertical= 1;
			else if (strcmp(next,"OFF")==0)
				temp->Print_Vertical= 0;
			else {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"OF")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, a filename is expected\n",argv[i]);
				exit(1);
			}
			strcpy(temp->Filename,argv[i+1]);
			name_set= 1;
			i+=2;
		}	
		else if (strcmp(inupper,"OP")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"ON")==0)
				temp->Prosite_Style= 1;
			else if (strcmp(next,"OFF")==0)
				temp->Prosite_Style= 0;
			else {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}	
		else if ((strcmp(inupper,"ON")==0) || (strcmp(inupper,"OA")==0)) {
			int par;
			if (!isnext) {
				fprintf (stderr,"After option %s, a non-negative integer value is expected\n",argv[i]);
				exit(1);
			}
			par=atoi(argv[i+1]);
			if (par<0) {
				fprintf (stderr,"After option %s, a non-negative integer value is expected\n",argv[i]);
				exit(1);
			}
			if (inupper[1]=='N')
				temp->Length_Hit_List=par;
			else
				temp->Number_Of_Occ_Output=par;
			i+=2;
		}
		else if ((strcmp(inupper,"PL")==0) || (strcmp(inupper,"PN")==0) || (strcmp(inupper,"PX")==0)) {
			int par;
			if (!isnext) {
            fprintf (stderr,"After option %s, a non-negative integer value is expected\n",argv[i]);
            exit(1);
         }
			par=atoi(argv[i+1]);
			if (par<0) {
				fprintf (stderr,"After option %s, a non-negative integer value is expected\n",argv[i]);
				exit(1);
			}
			switch (inupper[1]) {
				case 'L': temp->Max_Length_menu= par; break;
				case 'N': temp->Max_Num_Comp= max(1,par); break;
				case 'X': temp->Max_Gap= par; break;
			}
			i+=2;
		}
		else if (strcmp(inupper,"PF")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, a filename is expected\n",argv[i]);
				exit(1);
			}
			strcpy(temp->Restrictions_File,argv[i+1]);
			i+=2;
		}
		else if (strcmp(inupper,"PP")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, one of OFF, COMPLETE, or START is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"OFF")==0)
				temp->Restrictions_Input= off;
			else if (strcmp(next,"COMPLETE")==0)
				temp->Restrictions_Input= complete;
			else if (strcmp(next,"START")==0)
				temp->Restrictions_Input= start;
			else {
				fprintf (stderr,"After option %s, one of OFF, COMPLETE, or START is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"R")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"ON")==0)
				temp->Refinement= 1;
			else if (strcmp(next,"OFF")==0)
				temp->Refinement= 0;
			else {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"RG")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"ON")==0)
				temp->Refine_Generalise= 1;
			else if (strcmp(next,"OFF")==0)
				temp->Refine_Generalise= 0;
			else {
				fprintf (stderr,"After option %s, ON or OFF is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"S")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, one of INFO, MDL, DIST, TREE, or PPV is expected\n",argv[i]);
				exit(1);
			}
			if (strcmp(next,"INFO")==0) {
				temp->Tree_Input=0; temp->Diagnostic= 0; temp->Dist_Input= 0; temp->MDL_flag= 0;
			}
			else if (strcmp(next,"MDL")==0) {
				temp->Tree_Input=0; temp->Diagnostic= 0; temp->Dist_Input= 0; temp->MDL_flag= 1;
			}
			else if (strcmp(next,"DIST")==0) {
				temp->Tree_Input=0; temp->Diagnostic= 0; temp->Dist_Input= 1; temp->MDL_flag= 0;
			}
			else if (strcmp(next,"TREE")==0) {
				temp->Tree_Input=1; temp->Diagnostic= 0; temp->Dist_Input= 0; temp->MDL_flag= 0;
			}
			else if (strcmp(next,"PPV")==0) {
				temp->Tree_Input=0; temp->Diagnostic= 1; temp->Dist_Input= 0; temp->MDL_flag= 0;
			}
			else {
				fprintf (stderr,"After option %s, one of INFO, MDL, DIST, TREE, or PPV is expected\n",argv[i]);
				exit(1);
			}
			i+=2;
		}
		else if (strcmp(inupper,"SF")==0) {
			if (!isnext) {
				fprintf (stderr,"After option %s, a filename is expected\n",argv[i]);
				exit(1);
			}
			strcpy(temp->Tree_Filename,argv[i+1]);
			strcpy(temp->Dist_Filename,argv[i+1]);
			strcpy(temp->Swiss_Flat_File,argv[i+1]);
			i+=2;
		}
		else if ((strcmp(inupper,"Z0")==0) || (strcmp(inupper,"Z1")==0) ||
		         (strcmp(inupper,"Z2")==0) || (strcmp(inupper,"Z3")==0)) {
			float fpar;
			if (!isnext) {
				fprintf (stderr,"After option %s, a real is expected\n",argv[i]);
				exit(1);
			}
			fpar= atol(argv[i+1]);
			switch(inupper[1]) {
				case '0': temp->MDL_constant0= fpar;break;
				case '1': temp->MDL_constant1= fpar;break;
				case '2': temp->MDL_constant2= fpar;break;
				case '3': temp->MDL_constant3= fpar;break;
			}
			i+=2;
		}
		else {
			fprintf (stderr,"Unknown option %s\n",argv[i]);
			exit(1);
		}
	} /* i iterating over the arguments */
	temp->Max_Length=temp->Max_Length_menu+1;

#define MAX_MEMORY 3000000
	if (temp->Automatic) 
	{
	/* Automatically set the number of pattern components so that we do not get into
		memory trouble */
		int max_comp;

		max_comp= MAX_MEMORY/((Total_Seq_Length/8)*(temp->Max_Flex_Prod+max(temp->Nr_Symbol_Block,temp->Max_Gap+1)));
		temp->Max_Num_Comp=min(max_comp,min(temp->Max_Num_Comp,temp->Max_Length));
		printf ("max= %d\n",temp->Max_Num_Comp);
	}

	return temp;
}
	

t_Block_Options* Options_From_File(char* filename)
/******************************************************************/
{
	return NULL;
}	

void Options_Write_To_File(FILE *file,t_Block_Options *options)
/******************************************************************/
{
	char restrict[20];
	fprintf(file,"\n");
	fprintf(file,"                Pratt version 2.1 \n\n");
	fprintf(file,"        Analysing %d sequences from file %s\n\n",nrseqs,options->filename);

	fprintf(file,"PATTERN CONSERVATION:\n");
	fprintf(file,   "   CM: min Nr of Seqs to Match               %4d\n",options->Min_Nr_Seqs_Matching);
	fprintf(file,   "   C%%: min Percentage Seqs to Match         %5.1f\n",100.0*(float)options->Min_Nr_Seqs_Matching/(float)nrseqs);

	fprintf(file,   "\nPATTERN RESTRICTIONS :\n");
	fprintf(file,   "   PP: pos in seq [off,complete,start]   %8s\n",((options->Restrictions_Input==off) ? "off" : ((options->Restrictions_Input==start) ? "start" : "complete" )));
	if (options->Restrictions_Input!=off) 
		fprintf(file,"   PF: Restrictions filename %20s\n",(options->Restrictions_File));

	fprintf(file,   "   PL: max Pattern Length                    %4d\n",options->Max_Length_menu);
	fprintf(file,   "   PN: max Nr of Pattern Symbols             %4d\n",options->Max_Num_Comp);
	fprintf(file,   "   PX: max Nr of consecutive x's             %4d\n",options->Max_Gap);
	fprintf(file,   "   FN: max Nr of flexible spacers            %4d\n",options->Max_Num_Flex);
	if (options->Max_Num_Flex>0) {
		fprintf(file,"   FL: max Flexibility                       %4d\n",options->Max_Flex);
		fprintf(file,"   FP: max Flex.Product                      %4d\n",options->Max_Flex_Prod);
	}
	fprintf(file,   "   BI: Input Pattern Symbol File             %4s\n",(options->Input_Symbol_File ? "  on" : " off"));
	if (options->Input_Symbol_File) {
		fprintf(file,"   BF: Symbol Filename       %20s\n",options->Symbol_File);
	}
	fprintf(file,   "   BN: Nr of Pattern Symbols Initial Search  %4d\n",options->Nr_Symbol_Block);
	fprintf(file,"\nPATTERN SCORING:\n");
	fprintf(file,   "   S: Scoring [info,mdl,tree,dist,ppv]       %4s\n",(options->Tree_Input?"tree":(options->Dist_Input?"dist":(options->MDL_flag?"mdl":(options->Diagnostic?"ppv":"info")))));
	if (options->Tree_Input)
		fprintf(file,"   SF: Tree Filename         %20s\n",(options->Tree_Filename));
	if (options->Dist_Input)
		fprintf(file,"   SF: Distances Filename    %20s\n",(options->Dist_Filename));
	if (options->Diagnostic)
		fprintf(file,"   SF: Swiss-Prot Filename   %20s\n",(options->Swiss_Flat_File));

	/* Do not print MDL options on standard menu unless it is switched on 
      Want it in the menu summary when MDL_Pratt is used, but not for the
	   normal user - to avoid confusion                                    */
	if (options->MDL_flag) {
		fprintf(file,"MDL PARAMETERS (see ISMB96 paper):\n");
		fprintf(file,"   Z0:                                      %5.2f\n",options->MDL_constant0);
		fprintf(file,"   Z1:                                      %5.2f\n",options->MDL_constant1);
		fprintf(file,"   Z2:                                      %5.2f\n",options->MDL_constant2);
		fprintf(file,"   Z3:                                      %5.2f\n",options->MDL_constant3);
	}
	fprintf(file,"\nSEARCH PARAMETERS:\n");
	fprintf(file,   "   G: Pattern Graph from [seq,al,query]     %5s\n",(options->Alignment_flag ? "al" : (options->Query_Mode ? "query" : "seq")));
	if (options->Alignment_flag) {
		fprintf(file,"   GF: Alignment Filename    %20s\n",(options->Al_Filename));
	}
	else if (options->Query_Mode) {
		fprintf(file,"   GF: Query Sequence Filename %18s\n",options->Query_File_Name);
	}
	fprintf(file,   "   E: Search Greediness                      %4d\n",options->Quotient);
/*
	fprintf(file,   "   U: Pull-Out Parameter                     %4d\n",options->Quotient_Inc);
*/
	fprintf(file,   "   R: Pattern Refinement                     %4s\n",(options->Refinement ? "  on" : " off"));
	if (options->Refinement)
		fprintf(file,"   RG: Generalise ambiguous symbols          %4s\n",(options->Refine_Generalise ? "  on" : " off"));
	
	fprintf(file,"\n");
	fprintf(file,"OUTPUT:\n");
	fprintf(file,   "   OF: Output Filename       %20s\n",options->Filename);
	fprintf(file,   "   OP: PROSITE Pattern Format                %4s\n",(options->Prosite_Style ? "on" : "off"));
	fprintf(file,   "   ON: max number patterns                   %4d\n",options->Length_Hit_List);
	fprintf(file,   "   OA: max number Alignments                 %4d\n",options->Number_Of_Occ_Output);
	fprintf(file,   "   M: Print Patterns in sequences            %4s\n",(options->Print_Seqs_With_Motifs ? "on" : "off")); 
	if (options->Print_Seqs_With_Motifs)  {
		fprintf(file,"   MR: ratio for printing                    %4d\n",options->Print_Ratio);
		fprintf(file,"   MV: print vertically                      %4s\n",(options->Print_Vertical ? "on" : "off"));
	}
	fprintf(file,"\n");

}


int Menu_Summarise_Options(t_Block_Options *options)
/******************************************************************/
{
	printf ("\n");
	printf ("     PARTIAL SUMMARY OF THE SEARCH PARAMETERS\n");
	printf ("     ----------------------------------------\n");
	printf ("\n");
	printf ("You have chosen to search for patterns matching at least %d\n",options->Min_Nr_Seqs_Matching);
	printf ("of the %d sequences given.\n",nrseqs);
	printf ("\n");
	printf ("Pattern constraints:\n");
	printf ("  pattern length (length of A-x(4)-B is 6)      <= %2d \n",options->Max_Length-1);
	printf ("  number of pattern symbols                     <= %2d\n",options->Max_Num_Comp);
	printf ("    A and [DE] are symbols of A-x(2,4)-[DE]\n");
	printf ("  wildcard length (length of x(1,3) is 3)       <= %2d \n",options->Max_Gap);
	printf ("  number of flexible wildcards                  <= %2d \n",options->Max_Num_Flex);
	if (options->Max_Num_Flex>0) {
		printf ("  wildcard flexibility (flex. of x(1,3) is 2)   <= %2d \n",options->Max_Flex);
		printf ("  flexibility of product                        <= %2d \n",options->Max_Flex_Prod);
		printf ("    flex of prod of A-x(1,3)-B-x(3,4)-C is \n");
		printf ("    (3-1+1)*(4-3+1)=3*2=6\n");
	}
	printf ("\n");
	printf ("  pattern symbols during initial pattern search are restricted \n");
	if (options->Input_Symbol_File)
		printf ("  to the %d  first entries in the file %s\n", options->Nr_Symbol_Block,options->Symbol_File);
	else
		printf ("  to the %d  first entries in the default symbol set\n", options->Nr_Symbol_Block);

	printf ("\n");
	printf ("The degree of greediness of the initial search is %d\n",options->Quotient);
	printf ("(0 gives exhaustive search and bigger numbers give more greedy search).\n");
	printf ("\n");
	printf ("Maximum %d patterns will be reported, and for the %d best, the matching\n",options->Length_Hit_List,min(options->Number_Of_Occ_Output,options->Length_Hit_List));
	printf ("sequence segments will be shown aligned the pattern.\n",options->Number_Of_Occ_Output);
	printf ("\n");
	printf ("The output will be written to file %s\n",options->Filename);
	printf ("\n");
}

int OK_Parameter_Setting(t_Block_Options *options)
/******************************************************************/
{
	if ((options->Tree_Input) && (options->Dist_Input || options->MDL_flag || options->Diagnostic))
		return 0;
	if ((options->Dist_Input) && (options->Tree_Input || options->MDL_flag || options->Diagnostic))
		return 0;
	if ((options->MDL_flag) && (options->Tree_Input || options->Dist_Input || options->Diagnostic))
		return 0;
	if ((options->Diagnostic) && (options->Tree_Input || options->Dist_Input || options->MDL_flag))
		return 0;

	return 1;
}






/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
