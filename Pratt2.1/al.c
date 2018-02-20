/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>

#include "sequence.h"
#include "al.h"
#include "menu.h"

extern t_Block_Options *B_Options;

#define min(a,b) ((a)<(b)?(a):(b))


int PAM[26][26];

static int num_aa_line(char* line)
/****************************************************************************/
{
	int nr,pos;

	pos=16;
	nr= 0;
	while (isupper(line[pos]) || line[pos]=='-') {
		nr++;pos++;
	}
	return nr;
}

static t_alignment* Al_Get_Size(char* filename)
/****************************************************************************/
{
	FILE* file;
	int nr,l;
	char buffer[100];
	int finished;
	int n;
	t_alignment* al;

	if ((file=fopen(filename,"r"))==NULL) {
		printf ("File %s could not be opened\n", filename);
		exit(1);
	}

	al= (t_alignment*)calloc(1,sizeof(t_alignment));

	n=nr=l=0;

	if (fgets(buffer,100,file)==NULL) {
		printf ("Problems reading file %s\n",filename);
		exit (1);
	}

	if ((strstr(buffer,"CLUSTAL W"))==NULL) {
		printf ("First line of %s does not contain CLUSTAL W\n",filename);
		printf ("File needs to be in Clustal W format\n");
		exit (1) ;
	}

	fgets(buffer,100,file);
	while ((int)strlen(buffer)<5)
		fgets(buffer,100,file);
	
	finished= 0;
	while (!finished)
	{
		nr= 0;
		l+=num_aa_line(buffer);
		while (buffer[0]!=' ') {
			nr++;
			fgets(buffer,100,file);
		}

		if (n==0)
			n= nr;
		else
			assert (nr==n);

		if ((fgets(buffer,100,file))==NULL) 
			finished= 1;

		if ((fgets(buffer,100,file))==NULL) 
			finished= 1;
	}
	fclose(file);

	printf ("n=%d, l=%d\n",n,l);

	al->nr_seq= nr;
	al->nr_pos= l;

	return al;
}

static int Al_Get_Sequences(t_alignment* al,char* filename)
/****************************************************************************/
{
	FILE* file;
	int bulk;
	char buffer[100];
	int finished;
	int n;
	int i,j;
	int line= 0;

	if ((file=fopen(filename,"r"))==NULL) {
		printf ("File %s could not be opened\n");
		return 0;
	}

	al->positions= calloc(al->nr_seq,sizeof(char*));
	al->names= calloc(al->nr_seq,sizeof(char*));
	assert (al->positions && al->names);
	al->seq_index= NULL;

	for (i=0;i<al->nr_seq;i++) {
		al->positions[i]= calloc(al->nr_pos+1,sizeof(char));
		al->names[i]= calloc(100,sizeof(char));
		assert (al->positions[i] && al->names[i]);
	}

	fgets(buffer,100,file);
	line++;

	for (bulk=0;bulk< (al->nr_pos/60)+ ((al->nr_pos%60)?1:0) ;bulk++)
	{
		fgets(buffer,100,file);
		line++;
		while ((buffer[0]==' ') || (int)strlen(buffer)<5) {
			fgets(buffer,100,file);
			line++;
		}

		for (i=0;i<al->nr_seq;i++) {
			if (bulk==0)
				sscanf (buffer,"%s",al->names[i]);
			else
				if (strstr(buffer,al->names[i])!=buffer) {
					printf ("Inconsistency at line %d in %s: %s not like %s\n",
                        line,filename,al->names[i],buffer);
					exit (1);
				}
			j= 16;
			while (isupper(buffer[j])||buffer[j]=='-') {
				al->positions[i][bulk*60+(j-16)]= buffer[j];
				j++;
			}
			al->positions[i][bulk*60+(j-16)]='\0';
			fgets(buffer,100,file);
			line++;
		}
	}
	
	for (i=0;i<al->nr_seq;i++)
		printf ("%2d: %s\n",i,al->names[i]);

	for (i=0;i<al->nr_seq;i++) 
		assert (strlen(al->positions[i])==al->nr_pos);

	fclose(file);

	return 1;
}


t_alignment* Alignment_From_File (char* filename)
/****************************************************************************/
{
	t_alignment* al;

	printf ("The alignment in file %s is used to guide the search.\n\n",filename);

	if ((al= Al_Get_Size(filename))==NULL)
		return NULL;
	if (!(Al_Get_Sequences(al,filename)) )
		printf ("feil\n");

	return (al);	
}

t_alignment* Al_From_Sequences(t_sequence **s,int nr)
/****************************************************************************/
{
	t_alignment *al;
	int i;
	char *pt;
	int *ipt,*ept;

	al= malloc(sizeof(t_alignment));
	assert (al);
	al->nr_seq= 1;
	for (al->nr_pos=i=0;i<nr;i++)
		al->nr_pos+= s[i]->length;
	al->positions= malloc(sizeof(char*));
	al->positions[0]= malloc((al->nr_pos+10)*sizeof(char));
	assert (al->positions && al->positions[0]);
	al->seq_index=malloc(al->nr_pos*sizeof(int));
	assert (al->seq_index);
	al->names= malloc(sizeof(char*));
	al->names[0]= malloc(50*sizeof(char));
	assert (al->names && al->names[0]);
	if (nr==1)
		strcpy(al->names[0],s[0]->filename);
	else
		sprintf (al->names[0],"comb(%3d)",nr);
	for (pt=al->positions[0],i=0;i<nr;i++) {
		memcpy(pt,s[i]->seq,s[i]->length*sizeof(char));
		pt+= s[i]->length;
	}
	ipt= al->seq_index;
	for (i=0;i<nr;i++) {
		ept=ipt+s[i]->length;
		while (ipt<ept)
			*ipt++= i;
		/* definetely slower!!!!
		while (ipt<ept)
			*ipt++= 0;
		*/
	}
	return al;
}

static int Longer_Seq(const void *s1, const void *s2)
{
	return ((*(t_sequence**)s1)->length-(*(t_sequence**)s2)->length);
}

t_alignment* Alignment_From_Shortest_Sequence (int nr, t_sequence** s,int min_match)
/****************************************************************************/
{
	int i,l,index;
	t_alignment* al;
	t_sequence **ss;

	ss= malloc(sizeof(t_sequence*)*nr);
	assert (ss);
	memcpy(ss,s,nr*sizeof(t_sequence*));

	qsort((void*)ss,(size_t)nr,sizeof(t_sequence*),Longer_Seq);

	if (B_Options->Query_Mode)
		printf ("The special query sequence is used to guide the search.\n\n");
	if (nr==min_match)
		printf ("The shortest sequences is used to guide the search.\n\n");
	else if (nr>1)
		printf ("The %d shortest sequences are used to guide the search.\n\n",nr-min_match+1);

	/*
	printf ("The following %d sequences are used to guide the search:\n",nr-min_match+1);
	for (i=0;i<nr-min_match+1;i++)
		printf ("%s\n",ss[i]->filename);
	*/

	if (nr==1)
		al= Al_From_Sequences(ss,1);
	else
		al= Al_From_Sequences(ss,nr-min_match+1);

	assert (al);

	free (ss);
	return al;	
}

void Alignment_Print (t_alignment* al)
/****************************************************************************/
{
	int bulk;
	int i,j,k;
	int index;

	printf ("ALIGNMENT %d positions, %d sequences:\n",al->nr_pos,al->nr_seq);
	for (bulk=0;bulk< (al->nr_pos/60)+ ((al->nr_pos%60)?1:0) ;bulk++)
	{
		for (i=0;i<al->nr_seq;i++) {
			index= i;
			printf ("%10s ",al->names[index]);
			for (j=bulk*60;j<min((bulk+1)*60,al->nr_pos);j++)
				printf ("%c",al->positions[index][j]);
			printf ("\n");
		}
		printf ("\n");
	}
	printf ("\n\n");
}

void Alignment_Print_Seqs_Fasta(t_alignment* al,char *filename)
/****************************************************************************/
{
	FILE *fil;
	int i,j,k;
	int index;

	if ((fil= fopen(filename,"w"))==NULL)
	{
		fprintf (stderr,"Error opening %s in Alignment_Print_Seqs_Fasta\n",filename);
		exit(1);
	}

	for (i=0;i<al->nr_seq;i++) {
		fprintf (fil,"> %s\n",al->names[i]);
		k= 0;
		for (j=0;j<al->nr_pos;j++) {
			if (isupper(al->positions[i][j])) {
				fprintf (fil,"%c",al->positions[i][j]);
				k++;
				if ((k+1)%60==0)
					fprintf (fil,"\n");
			}
		}
		if ((k+1)%60!=0)
			fprintf (fil,"\n");
	}
	fclose(fil);
}
/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
