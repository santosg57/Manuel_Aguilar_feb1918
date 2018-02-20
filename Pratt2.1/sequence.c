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
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define round(a)  (floor((a)+0.5)) 

#define SEQ_MAIN
#include "sequence.h"

#include "menu.h"

#define MAX_SYMBOL 128

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))



static int Symbol_Initialized= 0;
static unsigned char DNA_Symbol[MAX_SYMBOL];
static unsigned char RNA_Symbol[MAX_SYMBOL];
static unsigned char PROTEIN_Symbol[MAX_SYMBOL];
extern FILE *mat;

double Prob_Random_Match;
int Min_Pos_Motif;
int Max_Ratio;

extern t_Block_Options* B_Options;

static void Initialize_Symbol_Table(char* alpha, unsigned char* table)
/**************************************************************************/
{
	int i;
	for (i=0;i<MAX_SYMBOL;i++)
		table[i]= 0;
	for (i=0;i<(int)strlen(alpha);i++)
		table[(int)alpha[i]]= 1;
}
		
static void Initialize_Symbol_Tables()
/**************************************************************************/
{
	if (!Symbol_Initialized)
	{	
		Symbol_Initialized= 1;
		Initialize_Symbol_Table(DNA_alphabet,DNA_Symbol);
		Initialize_Symbol_Table(RNA_alphabet,RNA_Symbol);
		Initialize_Symbol_Table(PROTEIN_alphabet,PROTEIN_Symbol);
	}
}



void seq_error_message (int error_nr,char* fil, int linje)
/**************************************************************************/
{
	printf ("Error in %s (line %d):",fil,linje);
	switch (error_nr) {
		case STORAGE_ERROR: 
			printf ("Storage allocation error\n");
			break;
		case FILE_OPEN_ERROR:
			printf ("File could not be opened\n");
			break;
		case ROTTEN_ERROR:
			printf ("There is something rotten in the state of Denmark!\n");
			break;
		default: 
			printf ("Unknown error\n");
			break;
	}
	exit (1);
}
	
/* LOCAL FUNCTIONS :*/
static void seq_expand(t_sequence *s)
/**************************************************************************/
{
	s->max_length*=2;
	if (!(s->seq=(char*)realloc(s->seq,sizeof(char)*s->max_length)))
		seq_error_message(1,__FILE__,__LINE__);
}

/* EXPORTED FUNCTIONS */

int Sequence_Initialize_Translation()
/**************************************************************************/
{
   FILE* fil;
   int i,j;
   char st[2];

   printf ("Sequence_Initialize_Translation\n");
   if ((fil= fopen("trans","r"))==NULL)
   {
		printf ("Could not open required file trans\n");
		seq_error_message (FILE_OPEN_ERROR,__FILE__,__LINE__);
   }

   for (i=0;i<128;i++)
      Trans_Table[i]= (char) 0;

   for (i=0;i<20;i++)
   {
      fscanf (fil,"%1s %d",st,&j);
      printf (":%c .%d\n",st[0],j); 
      if ( j<0 || j>=20 )
         printf ("error in line %d in file trans\n",i);
      Trans_Table[(int)st[0]]=(char)1+j;
   }
	fscanf (fil,"%lf",&Prob_Random_Match);
	printf ("Prob_Random_Match=%g\n", Prob_Random_Match);
	fscanf (fil,"%d",&Min_Pos_Motif);
	printf ("Min_Pos_Motif=%d\n", Min_Pos_Motif);
	fscanf (fil,"%d",&Max_Ratio);
	printf ("Max_Ratio=%d\n",Max_Ratio);

   fclose(fil);
  
   return 1;
}

	
t_sequence *seq_init_sequence(int nr,t_seq_type type,char* name)
/**************************************************************************/
{
	t_sequence *temp;

	temp= (t_sequence*)calloc(1,sizeof(t_sequence));
	temp->filename= (char*)calloc(strlen(name)+1,sizeof(char));
	strcpy(temp->filename,name);
	temp->length= 0;
	temp->max_length= nr;
	temp->seq=(char*)calloc(nr,sizeof(char));
	temp->type= type;
	temp->pat_position= (-1);
	temp->restricted= 0;

	return temp;
}

void seq_add_res(t_sequence *s,char t)
/**************************************************************************/
{
	if (s->length==s->max_length)
		seq_expand(s);
	s->seq[s->length]=t;
	s->length++;
}

static FILE* open_gcg_file_and_read_preamble(char* filename)
/**************************************************************************/
{
	FILE* fil;

	if (fil=fopen(filename,"r"))
	{
		int s;
		char c;

		s=0;
		while ( (s!=2) && 
				  ((char)(c=(char)fgetc(fil))!=(char)EOF) )
		{
			switch (s) 
			{
				case 0: 
					if (c=='.') 
						s= 1;
					break;
				case 1:
					s= (c=='.' ? 2 : 0);
					break;
				case 2:
				default:
					break;
			}
		} /* while */
		if (s==2)
			return (fil);
	} /* if */
	else
	{
		printf ("File %s could not be opened\n");
	}
	fclose (fil);
	return (NULL);
}

static int GCGfile_type(FILE* fil,unsigned char table[])
/**************************************************************************/
{
	char c;
	while ((char)(c=(char)fgetc(fil))!=(char)EOF)
	{
		if (isprint(c) && (c!=' ') && (!isdigit(c)) && isupper(c) && (table[(int)c]==0) )
		{
/*			printf ("->%c<- ",c); */
			return 0;
		}
	}
	return 1;
}

t_seq_type get_GCGfile_sequence_type(char* filename)
/**************************************************************************/
{
	FILE* fil;
	
	return protein;

	if ((fil=open_gcg_file_and_read_preamble(filename))!=NULL)
	{
		if (GCGfile_type(fil,DNA_Symbol))
			{fclose (fil);return dna;}
		else if (GCGfile_type(fil,RNA_Symbol))
			{fclose (fil);return rna;}
		else if (GCGfile_type(fil,PROTEIN_Symbol))
			{fclose (fil);return protein;}
		else
			{fclose (fil);return unknown;}
	}
	return unknown;
}


t_sequence *seq_from_GCG_file (char* filename)
/**************************************************************************/
{
	FILE* fil;
	t_sequence *seq;
	t_seq_type type;

	Initialize_Symbol_Tables();

	printf ("%s\n",filename);

	if ( ((type=get_GCGfile_sequence_type(filename))!=unknown)  &&
	     (fil=open_gcg_file_and_read_preamble(filename)) )
	{
		char c;
		seq=seq_init_sequence(400,type,filename);

		while ((char)(c=(char)fgetc(fil))!=(char)EOF)
		{
			if (isprint(c) && !isdigit(c) && c!=' ' && isupper(c))
				seq_add_res(seq,c);
		}
		fclose (fil);
		return (seq);
	}
	else
	{
		printf ("File %s could not be opened, or did not contain valid sequence\n",filename);
		return (NULL);
	}
}

t_sequence* one_seq_from_Fasta_file(char *filename,t_seq_type type)
/**************************************************************************/
{
	FILE* fil;
	t_sequence *seq;
	char buf[100];
	char name[20];
	int line;
	int max_nr;
	int i;
	char *ok;

   seq= NULL;

	Initialize_Symbol_Tables();

	fil=fopen(filename,"r");
	if (fil==NULL) {
		printf ("Could not open file %s\n",filename);
		return NULL;
	}
	
	line= 0;
	ok= buf;
	ok= fgets(buf,100,fil);
	while (ok!=NULL)
	{
		while (buf[0]!='>' && ok ) {
			ok= fgets(buf,100,fil);
			line++;
		}
		if (!ok)
			continue;
		sscanf(buf+1,"%s",name);
		if (seq!=NULL)
		{
			fprintf (stderr,"read query file %s: contains more than one sequence\n",filename);
			fprintf (stderr,"using first\n");
			return seq;
		}
		seq=seq_init_sequence(400,type,name);
		while ( (ok=fgets(buf,100,fil)) && (buf[0]!='>' ) )
		{
			line++;
			for (i=0;i<(int)strlen(buf);i++)
				if (isupper(buf[i]) || islower(buf[i]) )
				{
					if (islower(buf[i]))
						buf[i]= toupper(buf[i]);
					seq_add_res(seq,buf[i]);
				}
		}
	}
	return seq;
}


t_sequence** seqs_from_Fasta_file(char *filename,int *nr,t_seq_type type)
/**************************************************************************/
{
	FILE* fil;
	t_sequence **seqs;
	t_sequence *seq;
	char buf[100];
	char name[20];
	int line;
	int max_nr;
	int i;
	char *ok;

	Initialize_Symbol_Tables();


	fil=fopen(filename,"r");
	if (fil==NULL) {
		printf ("Could not open file %s\n",filename);
		return NULL;
	}
	max_nr= 10000;
	seqs= (t_sequence**)calloc(max_nr,sizeof(t_sequence*));
	
	line= 0;
	(*nr)= 0;
	ok= buf;
	ok= fgets(buf,100,fil);
	while (ok!=NULL)
	{
		while (buf[0]!='>' && ok ) {
			ok= fgets(buf,100,fil);
			line++;
		}
		if (!ok)
			continue;
		sscanf(buf+1,"%s",name);
		seq=seq_init_sequence(400,type,name);
		while ( (ok=fgets(buf,100,fil)) && (buf[0]!='>' ) )
		{
			line++;
			for (i=0;i<(int)strlen(buf);i++)
				if (isupper(buf[i]) || islower(buf[i]) )
				{
					if (islower(buf[i]))
						buf[i]= toupper(buf[i]);
					seq_add_res(seq,buf[i]);
				}
		}
		seqs[(*nr)]= seq;
		(*nr)++;
	}
	return seqs;
}

t_sequence** seqs_from_SwissProt_file(char *filename,int *nr,t_seq_type type)
/**************************************************************************/
{
	FILE* fil;
	t_sequence **seqs;
	t_sequence *seq;
	char buf[100];
	char name[20];
	int line;
	int max_nr;
	int i;
	char *ok;

	Initialize_Symbol_Tables();


	fil=fopen(filename,"r");
	if (fil==NULL) {
		printf ("Could not open file %s\n",filename);
		return NULL;
	}
	max_nr= 10000;
	seqs= (t_sequence**)calloc(max_nr,sizeof(t_sequence*));
	assert (seqs);
	
	line= 0;
	(*nr)= 0;
	ok= filename;
	while (ok!=NULL)
	{
		while ((ok=fgets(buf,100,fil)) && (buf[0]!='I' || buf[1]!='D'))
			line++;
		if (!ok)
			continue;
		sscanf(buf+5,"%s",name);
		seq=seq_init_sequence(400,type,name);
		while ( (ok=fgets(buf,100,fil)) && (buf[0]!='S' || buf[1]!='Q'))
			line++;
		if (!ok)
			continue;
		while ( (ok=fgets(buf,100,fil)) && (buf[0]!='/' || buf[1]!='/'))
		{
			line++;
			for (i=0;i<(int)strlen(buf);i++)
				if (isupper(buf[i]))
					seq_add_res(seq,buf[i]);
		}
		if (!ok) {
			fprintf (stderr,"sequence %d not terminated\n",(*nr)+1);
			return NULL;
		}
		seqs[(*nr)]= seq;
		(*nr)++;
	}
	return seqs;
}


t_sequence** seqs_from_SRS_file(char *filename,int *nr,t_seq_type type)
/**************************************************************************/
{
	FILE* fil;
	t_sequence **seqs;
	t_sequence *seq;
	char buf[100];
	char name[20];
	int line;
	int max_nr;
	int i;
	char *ok;

	Initialize_Symbol_Tables();


	fil=fopen(filename,"r");
	if (fil==NULL) {
		printf ("Could not open file %s\n",filename);
		return NULL;
	}
	max_nr= 10000;
	seqs= (t_sequence**)calloc(max_nr,sizeof(t_sequence*));
	assert (seqs);
	
	line= 0;
	(*nr)= 0;
	ok= filename;
	while (ok!=NULL)
	{
		while ((ok=fgets(buf,100,fil)) && (buf[0]!='I' || buf[1]!='D'))
			line++;
		if (!ok)
			continue;
		sscanf(buf+5,"%s",name);
		seq=seq_init_sequence(400,type,name);
		while ( (ok=fgets(buf,100,fil)) && (buf[0]!='S' || buf[1]!='Q'))
			line++;
		if (!ok)
			continue;
		while ( (ok=fgets(buf,100,fil)) && (strlen(buf)>1))
		{
			line++;
			if (buf[0]=='>' || (buf[0]==' '&& strstr(buf,"Length")))
				continue;
			for (i=0;i<(int)strlen(buf);i++)
				if (isupper(buf[i]))
					seq_add_res(seq,buf[i]);
		}
		if (!ok) {
			printf ("sequence %d not terminated\n",(*nr)+1);
			return NULL;
		}
		seqs[(*nr)]= seq;
		(*nr)++;
	}
	return seqs;
}



void seq_print_seq(t_sequence* seq,int nr)
/**************************************************************************/
{
	int i;

	printf ("sequence %2d: %s (",seq->index,seq->filename);
	switch(seq->type) {
		case dna: printf ("DNA"); break;
		case rna: printf ("RNA"); break;
		case protein: printf ("protein"); break;
		case unknown: printf ("unknown"); break;
	}
	printf (",%d)",seq->length);
	for (i=1;i<=seq->length;i++)
	{
		if (((i-1)%nr)==0) 
			printf("\n%5d ",i);
		if ((((i-1)%nr)%10)==0)
			printf (" ");
		printf("%c",seq_element(seq,i-1));
	}
	printf("\n\n");
}

void seq_print_part_file(FILE* file,t_sequence* seq, int i1, int i2)
/**************************************************************************/
{
	int i;

#ifdef DEBUG
	assert (i2>=i1);
	assert (i1>=0);
#endif
	
	if (i2>=seq->length)
		i2= seq->length-1;

	fprintf (file,"%20s:%6d-%6d: ",seq->filename,i1+1,i2+1);
	for (i=i1-5;i<i1;i++)
		fprintf (file,"%c",(i>=0 ? (tolower(seq_element(seq,i))) : ' ') );
	fprintf (file," ");
	for (i=i1;i<=i2;i++)
		fprintf(file,"%c",seq_element(seq,i));
	fprintf (file," ");
	for (i=i2+1;i<=i2+5;i++)
		fprintf (file,"%c",(i<seq->length ? (tolower(seq_element(seq,i))) : ' ') );
	fprintf(file,"\n");
}

void seq_print_part(t_sequence* seq, int i1, int i2)
/**************************************************************************/
{
	seq_print_part_file(stdout,seq,i1,i2);
}

t_seq_freq* seq_freq_from_seq(t_sequence* seq)
/**************************************************************************/
{
	t_seq_freq* temp;
	char* alpha;
	int i,j;

	temp= (t_seq_freq*)calloc(1,sizeof(t_seq_freq));

	switch (seq->type) {
	case dna:
		temp->nr= 4;
		alpha= DNA_alphabet;
		break;
	case rna:
		temp->nr= 4;
		alpha= RNA_alphabet;
		break;
	case protein:
		temp->nr=20;
		alpha=PROTEIN_alphabet;
		break;
	case unknown:
		seq_error_message(3,__FILE__,__LINE__);
		break;
	}
	temp->base=(char*)calloc(temp->nr,sizeof(char));
	temp->freq=(int*)calloc(temp->nr,sizeof(int));
	for (i=0;i<temp->nr;i++)
	{
		temp->freq[i]=0;
		temp->base[i]= alpha[i];
		for (j=0;j<seq->length;j++)
			if (seq->seq[j]==temp->base[i])
				temp->freq[i]++;
	}	
	for (i=0;i<temp->nr;i++)
	{
		temp->ratio[(int)temp->base[i]]=
			(float)(temp->freq[i])/(float)(seq->length);
	}

	temp->exp= 0.0;
	for (i=0;i<temp->nr;i++)
	{
		float t;
		t= (float)(temp->freq[i])/(float)(seq->length);
		temp->exp+= t*t;
	}
	printf ("exp:%f\n",temp->exp);

	return temp;
}

void seq_freq_print(t_seq_freq* freq)
/**************************************************************************/
{
	int i;
	int total;
	for (i=total=0;i<freq->nr;i++)
	{
		printf ("%c:     %5d\n",freq->base[i],freq->freq[i]);
		total+=freq->freq[i];
	}
	printf ("------------\ntotal: %5d\n",total);
}

/*
#define MAX_LEN_SIM 10000
static int sim[MAX_LEN_SIM+1];
float seq_dist (char* s1, int l1, char *s2, int l2)
{
	int i,j;
	int prev, eq, temp;

#ifdef DEBUG
	assert (l1<MAX_LEN_SIM);
#endif

	for (i=0;i<=l1;i++)
		sim[i]= 0;

	for (i=0;i<=l2;i++) {
		prev= 0;
		for (j=0;j<l1;j++) {
			eq= (s1[j]==s2[i] ? 1 : 0 );
			temp= sim[j+1];
			sim[j+1]= max ( prev + eq, max (sim [ j ], sim [j+1]));
			prev= temp;
		}
	}
	return ((float)sim[l1])/(float)(min(l1,l2));
}
*/

void Spec_Print_Sequence_To_File(FILE** file,t_sequence* seq)
/**************************************************************************/
{
	int i;

	fprintf (*file,"%s  ",seq->filename);
	for (i=0;i<seq->length;i++)
	{
		fprintf(*file,"%c",seq_element(seq,i));
	}
	fprintf (*file,"\n");
}

void Print_Sequence_To_File(FILE** file,t_sequence* seq)
/**************************************************************************/
{
	int i;

	fprintf (*file,"> %s\n",seq->filename);
	for (i=0;i<seq->length;i++)
	{
		fprintf(*file,"%c",seq_element(seq,i));
		/*
		if (i>0 && (i%60)==0)
		*/
		if (((i+1)%60)==0)
			fprintf (*file,"\n");
	}
	fprintf (*file,"\n");
}



void Print_Seqs_To_Fasta_File(t_sequence** sequences,int nr,char *filename)
/**************************************************************************/
{
	FILE *file;
	int i;

	file= fopen(filename,"w");

	if (file==NULL) {
		printf ("Cannot write Fasta format to %s\n",filename);
		return;
	}

	for (i=0;i<nr;i++)
		Print_Sequence_To_File(&file,sequences[i]);

	fclose (file);
}


void Seqs_Prepare_For_Output(t_sequence** sequences, int nr)
/**************************************************************************/
{
	int i,l;
	t_sequence *s;

	for (i=0;i<nr;i++)
	{
		s= sequences[i];
		l= s->length/(B_Options->Print_Ratio)+1;
		s->out= calloc(l+1,sizeof(char));
		memset(s->out,(int)'.',l);
	}
}

int Seq_Add_Symbol_To_Output(t_sequence* s, int pos, char symbol)
/**************************************************************************/
{
	if (s->out[pos/(B_Options->Print_Ratio)]=='.') {
		s->out[pos/(B_Options->Print_Ratio)]= symbol;
		return 1;
	}
	else {
		return (s->out[pos/(B_Options->Print_Ratio)]==symbol);
	}
}

#define LINE 60

static int Seqs_Print_Out_Vertical(FILE **file,t_sequence** sequences, int nr)
/**************************************************************************/
{
	int i,j,l,k, doit;
	t_sequence *s;
	char line[LINE+20];
	char blanks[LINE+2];
	char name[20];
	char *lptr;
	int remaining;
	int first;
	int more;
	int position;

	for (i=0;i<LINE;i++)
		blanks[i]=' ';
	blanks[LINE]='\0';

	fprintf (*file,"\n");
	for (i=0;i<LINE;i++)
		fprintf (*file,"-");
	fprintf (*file,"\n");

	remaining= nr;
	first= 0;
	while (remaining>0) {
		fprintf (*file,"#NEWPAGE\n");
		fprintf (*file,"matches in sequences %d to %d\n\n",first+1,first+1+min(remaining,LINE)-1);
		fprintf (*file,"seq   Sequences vertically\n");
		fprintf (*file,"pos.\n");
		/* PRINT NAMES */
		for (i=0;i<16;i++){
			strcpy(line,blanks);
  			lptr = line;
  			for (j=0;j<min(remaining,LINE);j++) {
    			s= sequences[first+j];
    			sprintf (name,"%16s",s->filename);
    			*lptr++ = name[i];
    		}
  			fprintf (*file,"      %s\n",line);
		}
		/* PRINT SEQUENCES */
		strcpy(line,blanks);
  		fprintf (*file,"      %s\n",line);
		position= 1;
		i=0; more= 1;
		while (more) {
			more= 0;
			strcpy(line,blanks);
  			lptr = line;
  			for (j=0;j<min(remaining,LINE);j++) {
    			s= sequences[first+j];
				if (i<(int)strlen(s->out)) {
					*lptr= s->out[i];
					more=1;
				}
				lptr++;
			}
			if (more) {
				if ((i%10)==0)
  					fprintf (*file,"%4d  %s\n",position,line);
				else
  					fprintf (*file,"      %s\n",line);
    		}
			i++;
			position+= B_Options->Print_Ratio;
		}
		fprintf (*file,"\n");
		remaining-= min(LINE,remaining);
		first+= LINE;
	}
	fprintf (*file,"\n");
	for (i=0;i<LINE;i++)
		fprintf (*file,"-");
	fprintf (*file,"\n");
}


int Seqs_Print_Out(FILE **file,t_sequence** sequences, int nr)
/**************************************************************************/
{
	int i,j,l;
	t_sequence *s;

	if (B_Options->Print_Vertical)
		Seqs_Print_Out_Vertical(file,sequences,nr);
	else
	{
		fprintf (*file,"\n-------------------------------------------------\n");
		for (i=0;i<nr;i++)
		{
			s= sequences[i];
			fprintf (*file,"%15s: ",s->filename);
			fprintf (*file,"%s\n",s->out);
		}
		fprintf (*file,"\n-------------------------------------------------\n");
	}
}


int Sequences_Read_Restrictions_File(t_sequence **ss,int nrs,char *filename)
/**************************************************************************/
{
	FILE *file;
	char buf[200];
	char name[100];
	char rest[100];
	int first[MAX_RESTRICT],last[MAX_RESTRICT];
	int line;
	int i,j;
	int nr;
	int nr_assigned;

	if ((file= fopen(filename,"r"))==NULL) 
		return 0;

	line= 0;
	while (fgets(buf,200,file)) {
		line++;
		if (buf[0]=='>') {
			if (sscanf(&(buf[1]),"%s %100c",name,rest)<2) {
				fprintf (stderr,"Error in line %d of restrictions file %s\n",line,filename);
				exit (1);
			}
			strcpy(buf,rest);
			nr= 0;
			while (strlen(buf)>1 && nr<MAX_RESTRICT) {
				memset(rest,0,100);
				if ((nr_assigned=sscanf(buf,"(%d,%d) %100c",&(first[nr]),&(last[nr]),rest))>=2) {
					if (last[nr]!=0 && last[nr]<first[nr]) {
						fprintf (stderr,"Error in line %d of restrictions file %s:",line,filename);
						fprintf (stderr,"allowed range (%d,%d) is empty\n",first[nr],last[nr]);
						exit (1);
					}
					if (first[nr]<0) {
						fprintf (stderr,"Error in line %d of restrictions file %s:",line,filename);
						fprintf (stderr,"negative numbers in range (%d,%d)\n",first[nr],last[nr]);
						exit (1);
					}	
					nr++;
					if (nr_assigned==3)
						strcpy(buf,rest);
					else
						memset(buf,0,100);
				}
				else 
					memset(buf,0,100);
			}
			if (strlen(buf)>1) { 
				fprintf(stderr,"**********************************\n");
				fprintf(stderr,"WARNING: not more than %d intervals are allowed -- check line %d of restrictions file %s\n",MAX_RESTRICT,line,filename);
				fprintf(stderr,"Thus only the first %d intervals given are used.\n",MAX_RESTRICT);
				fprintf(stderr,"**********************************\n");
				fprintf(stderr,"\n");
			}
			if (nr==0) {
				fprintf(stderr,"\nERROR in restrictions file\n");
				fprintf(stderr,  "--------------------------\n");
				fprintf(stderr,"Line %d in file %s with restrictions for sequence %s has wrong format.\n",line,filename);
				fprintf(stderr,"The format of restrictions is \n");
				fprintf(stderr,">[sename] ([lower],[upper])\n");
				fprintf(stderr,"i.e. lines starting with > followed by a sequence name and pairs of integers\n");
				fprintf(stderr,"enclosed in paranthesis. Several pairs (up to %d) can be given for each\n",MAX_RESTRICT);
				fprintf(stderr,"sequence, and for each sequence listed, at least one pair must be given.\n");
				fprintf(stderr,"\n");
				return 0;
			}
			i=0;
			while (i<nrs && (strcmp(name,ss[i]->filename)!=0))
				i++;
			if (i<nrs) {
				/* subtract one to get 'C' array indexes instead of residue numbers */
				for (j=0;j<nr;j++) {
					ss[i]->first[j]= max(0,min(first[j]-1,ss[i]->length-1));
					ss[i]->last[j]= (last[j]==0 ? ss[i]->length-1: min(last[j]-1,ss[i]->length-1));
				}
				ss[i]->restricted= nr;
			}
			else {
				fprintf(stderr,"Filename %s at line %d of restrictions file %s is not name of a sequence\n",name,line,filename);
				exit (1);
			}
		}
		else {
			fprintf(stderr,"\nERROR in restrictions file\n");
			fprintf(stderr,  "--------------------------\n");
			fprintf(stderr,"Line %d in file %s with restrictions for sequence %s has wrong format.\n",line,filename);
			fprintf(stderr,"The format of restrictions is \n");
			fprintf(stderr,">[sename] ([lower],[upper])\n");
			fprintf(stderr,"i.e. lines starting with > followed by a sequence name and pairs of integers\n");
			fprintf(stderr,"enclosed in paranthesis. Several pairs (up to %d) can be given for each\n",MAX_RESTRICT);
			fprintf(stderr,"sequence, and for each sequence listed, at least one pair must be given.\n");
			fprintf(stderr,"\n");
			return 0;
		}
	}
	return 1;
}

void Print_Restrictions(FILE* file,t_sequence **ss,int nrs)
/**************************************************************************/
{
	int i,j;

	fprintf(file,"Restrictions on pattern matches:\n");
	fprintf(file,"-------------------------------\n");
	for (i=0;i<nrs;i++) {
		if (ss[i]->restricted>0) {
			fprintf(file,"%14s:",ss[i]->filename);
			for (j=0;j<ss[i]->restricted;j++)
				fprintf(file,"(%d,%d)",ss[i]->first[j]+1,ss[i]->last[j]+1);
			fprintf(file,"\n");
		}
	}
	fprintf(file,"\n");
}

void Print_Seq_Lengths(FILE* file,t_sequence **ss,int nrs)
/**************************************************************************/
{
	int i;

	fprintf (file,"Sequence lengths:\n");
	for (i=0;i<nrs;i++)
		fprintf(file,"%14s %6d\n",ss[i]->filename,ss[i]->length);
}

/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
