/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "sequence.h"
#include "al.h"
#include "pattern.h"
#include "hit.h"
#include "block.h"
#include "menu.h"
#include "al.h"
#include "help.h"

#define BUFFER 20


#define MIN_ID 0.8
#define MAX_ID 0.5
#define MAX_NR_SEQS 1000

char version[]="2.1";

int nrseqs;
int virt_nrseqs;
t_sequence** seqs;
t_sequence* query_sequence;
int Total_Seq_Length;
time_t start_sec;
char start_time[32];


/*
FILE* stat;
*/
extern float **Features;
extern int Nr_Features;

void usage()
/**************************************************************************/
{
	printf ("usage: pratt <format> <filename> [options]\n");
	printf ("   where format is one of   swissprot, fasta, gcg\n");
	printf ("   and the file is on the given format\n");
	general_help();
	exit (1);
}


static int Initialize_From_File_List(char* filename)
/**************************************************************************/
{
	FILE* fil;
	char Line[80];
	char name[80];

	nrseqs= 0;
	seqs= (t_sequence**) calloc(MAX_NR_SEQS,sizeof(t_sequence*));

	if ((fil=fopen(filename,"r"))==NULL)
		return 0;

	while (fgets(Line,80,fil)!=NULL)
	{
		name[0]= '\0';
		if ( ((int)strlen(Line)!=0) && sscanf(Line,"%s",name) && ((int)strlen(name)>0))
		{
			if ( (seqs[nrseqs]= seq_from_GCG_file (name))==NULL)
				printf ("Error in reading file %s (%s,%d)\n",
						name,__FILE__,__LINE__);
			else 
			{
			/*
				if ((int)strlen(Line)>(int)strlen(name)+1) {
					sscanf (Line,"%s %d",name,&(seqs[nrseqs]->pat_position));
					seqs[nrseqs]->pat_position--;
				}
				else
					seqs[nrseqs]->pat_position= (-1);
			*/
				nrseqs++;
			}
		}
	}

	assert (nrseqs<MAX_NR_SEQS);
	fclose(fil);

	return 1;
}

void Print_Pratt_Info(FILE *file)
/**************************************************************************/
{
   fprintf (file,"\n");
   fprintf (file,"------------------------------------------------------------\n");
   fprintf (file,"                Pratt version %s, Febr. 1997\n",version);
   fprintf (file,"                 Written by Inge Jonassen, \n");
   fprintf (file,"                   University of Bergen\n");
   fprintf (file,"                           Norway\n");
   fprintf (file,"                   email: inge@ii.uib.no\n");
   fprintf (file,"                 For more information, see \n");
   fprintf (file,"            http://www.ii.uib.no/~inge/Pratt.html\n");
   fprintf (file,"------------------------------------------------------------\n");
   fprintf (file,"                        Please quote:\n");
   fprintf (file,"             I.Jonassen, J.F.Collins, D.G.Higgins.\n");
   fprintf (file,"             Protein Science 1995;4(8):1587-1595.\n");
   fprintf (file,"------------------------------------------------------------\n");
   fprintf (file,"\n\n");
}


main(int argc,char** argv)
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
{
	int i,j,k;
	int ok;
	char Pat_Filename[100];
	t_Block_Options *Options;
	char buf[20];
	FILE *out;
	char *local_time;

	Print_Pratt_Info(stdout);
	

	if (argc<3)
		usage();

	for (i=0;i<(int)strlen(argv[1]);i++)
		if (isupper(argv[1][i]))
			argv[1][i]= tolower(argv[1][i]);

	if ((strcmp(argv[1],"gcg"))==0) {
		if (!Initialize_From_File_List(argv[2])) {
			printf ("%s not a list of filenames\n",argv[2]);
			usage ();
		}
	}
	else if ((strcmp(argv[1],"srs"))==0) {
		if ((seqs=seqs_from_SRS_file(argv[2],&nrseqs,protein))==NULL) {
			printf ("%s is not a valid file\n",argv[2]);
			usage ();
		}
	}
	else if ((strcmp(argv[1],"swissprot"))==0) {
		if ((seqs=seqs_from_SwissProt_file(argv[2],&nrseqs,protein))==NULL) {
			printf ("%s is not a valid file\n",argv[2]);
			usage ();
		}
	}
	else if ((strcmp(argv[1],"fasta"))==0) {
		if ((seqs=seqs_from_Fasta_file(argv[2],&nrseqs,protein))==NULL) {
			printf ("%s is not a valid file\n",argv[2]);
			usage ();
		}
	}
	else {
		printf ("unknown file format\n");
		usage();
	}

/*
	if (argc>3) { 
		QUERY=1;
		query_sequence= one_seq_from_Fasta_file(argv[3],protein);
	}
*/

	for (Total_Seq_Length=i=0;i<nrseqs;i++) {
		seqs[i]->index= i+1;
		seqs[i]->restricted= 0;
		seqs[i]->first[0]= 0;
		seqs[i]->last[0]= seqs[i]->length-1;
		Total_Seq_Length+= seqs[i]->length;
	}	

	sprintf (Pat_Filename,"%s\0",argv[2]);

	if (nrseqs==0) {
		printf ("No sequences!!\n");
		exit (1);
	}

/*
	if (nrseqs==1) {
		printf ("Only one sequence, I dont bother!\n");
		exit (1);
	}
*/

/*
#define FASTA_OUT
*/
#ifdef FASTA_OUT
{
	char fasta_filename[100];
	sprintf (fasta_filename,"%s.fasta",argv[2]);
	/*
	sprintf (fasta_filename,"fasta/%s",argv[2]);
	*/
	Print_Seqs_To_Fasta_File(seqs,nrseqs,fasta_filename);
}
#endif


/*
	fprintf (stderr,"\n\n%s ",argv[2]);

	if (QUERY) {
		seqs[nrseqs]= query_sequence;
		nrseqs++;
	}
*/

	if (argc>3) {
		Options= Options_From_Argumentline(nrseqs,Pat_Filename,argc,argv);
		if (Options->Show_Menu) {
			Options= Menu_Set_Options(nrseqs,Pat_Filename,Options);
			Menu_Summarise_Options(Options);
			printf ("Are you happy with these settings (y/n)? ");
			scanf ("%s",buf);
			while  (buf[0]!='y' && buf[0]!='Y') {
				Options= Menu_Set_Options(nrseqs,Pat_Filename,Options);
				Menu_Summarise_Options(Options);
				printf ("Are you happy with these settings (y/n)? ");
				scanf ("%s",buf);
			}
		}
	}
	else	{
		Options= Menu_Set_Options(nrseqs,Pat_Filename,NULL);
		Menu_Summarise_Options(Options);
		printf ("Are you happy with these settings (y/n)? ");
		scanf ("%s",buf);
		while  (buf[0]!='y' && buf[0]!='Y') {
			Options= Menu_Set_Options(nrseqs,Pat_Filename,Options);
			Menu_Summarise_Options(Options);
			printf ("Are you happy with these settings (y/n)? ");
			scanf ("%s",buf);
		}
	}
	if (!OK_Parameter_Setting(Options)) {
		fprintf(stderr,"Illegal parameter settings: you cannot use multiple ranking systems simultaneously\n");
		fprintf(stderr,"Please rerun program with legal parameter sets\n");
		exit (1);
	}
	printf ("\n");
	printf ("Pattern search is starting\n");
	printf ("\n");
	
	start_sec= time(NULL);
	local_time = ctime(&start_sec);
	strcpy(start_time,local_time);

   if ((out= fopen(Options->Filename,"w"))==NULL)
   {
      printf ("Pratt cannot open file %s for writing \n",Options->Filename);
      printf ("Change write access or use another filename and run Pratt again\n");
      printf ("Ciao\n");
   }
	Print_Pratt_Info(out);
   Options_Write_To_File(out,Options);

	fprintf (out,"\n");
	fprintf (out,"\n");

	if (Options->Restrictions_Input) {
		if (!Sequences_Read_Restrictions_File(seqs,nrseqs,Options->Restrictions_File)) {
			printf ("Restrictions file %s could not be read or contained errors\n",Options->Restrictions_File);
			exit (1);
		}
	}

	if (Options->Restrictions_Input)  {
		fprintf (out,"\n");
		Print_Restrictions(out,seqs,nrseqs);
		fprintf (out,"\n");
	}

	if (Options->WWW) {
		fprintf (out,"\n");
		Print_Seq_Lengths(out,seqs,nrseqs);
		fprintf (out,"\n");
	}


	fprintf(out,"Pratt run started at %s\n",start_time);

   fclose(out);


	if (Options->Query_Mode)
		query_sequence= one_seq_from_Fasta_file(Options->Query_File_Name,protein);

	Find_Motif_Block(Options);

	exit(0);
}
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
