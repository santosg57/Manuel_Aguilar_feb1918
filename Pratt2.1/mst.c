/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
/*
Program to compute the minimum spanning tree of a complete
graph represented by the incidence matrix.
*/

#include<stdio.h>
#include "pam_dist.h"

#define MAX_NUM_DIST 500

#define inf 999

#define null -1

#define true 1
#define false 0



/*********************************************************/

static void print_1_int_array(int a[MAX_NUM_DIST], int imax) {

   int i,j;

   for (i=0; i<imax; i++) {
      printf("%6d ", a[i]);
   }
   printf("\n");
}

/*********************************************************/
  
static void print_1_float_array(float a[MAX_NUM_DIST], int imax) {

   int i,j;

   for (i=0; i<imax; i++) {
      printf("%6.3f ", a[i]);
   }
   printf("\n");
}

/*********************************************************/
  
static void print_2_float_array(float a[MAX_NUM_DIST][MAX_NUM_DIST], int imax, int jmax){

   int i,j;

   for (i=0; i<imax; i++) {
      for(j=0; j<jmax; j++) {
         printf("%6.3f ", a[i][j]);
      }
      printf("\n");
   }
}


/*********************************************************/

static int extract_min(int Q[MAX_NUM_DIST], float key[MAX_NUM_DIST], int imax){

   int i, istart, min_node ;

   istart = 0;
   while(Q[istart] == false) istart++;

   min_node = istart;

   for(i=istart+1; i<imax; i++){ 
      if((Q[i] == true) && (key[min_node] > key[i]))   
         min_node = i;
   }  

   return(min_node);
}


/*********************************************************/

float  mst_weight(float w[MAX_NUM_DIST][MAX_NUM_DIST], int size, int sub[MAX_NUM_DIST], char mst[])
{
   int Q[MAX_NUM_DIST], parent[MAX_NUM_DIST];
   float  key[MAX_NUM_DIST], minsum;

   int r;
   int i, j, u, v, QN, initial;

   char new_edge[20];
 
   QN = 0;
   for(i=0; i<size; i++) {
      Q[i] = sub[i];
      if(Q[i] == true) QN++;
      key[i] = inf;
      parent[i] = null;
   }

   r = 0;
   while((sub[r] == 0) && (r < size)) r++;

   key[r] = 0;
   mst[0] = '\0';
   minsum = 0;
  
   initial = true; 
   while( QN > 0){
      u = extract_min(Q, key, size);
      minsum += key[u];
      if(initial == false){
         sprintf(new_edge, "(%d,%d):%6.3f, ", parent[u], u, key[u]);
         strcat(mst, new_edge);
      }
     
/* tables for testing algorithm...
 
      printf("\n");
      printf("QN = %d\n", QN);
      printf("Q   :"); print_1_int_array(Q, size);
      printf("key :"); print_1_float_array(key, size);
      printf("par.:"); print_1_int_array(parent, size);
      printf("==>>  u = %d, minsum = %6.3f\n", u, minsum);

*/
      initial = false;
      Q[u] = false;
      QN--;

      for(v=0; v<size; v++){
         if((Q[v] == true)  && (w[u][v] < key[v])){
            parent[v] = u;
            key[v] = w[u][v];
         }
      }
   }     
/*
	printf ("mst tree:%s\n",mst);
*/
   return minsum;
}  
   

/*********************************************************/

int read_mat(float w[MAX_NUM_DIST][MAX_NUM_DIST],char *filename) {

   char text[200];
   int i, j, size;
	float temp;

   FILE *distfile;

   if ((distfile = fopen(filename, "r"))==NULL) {
		printf ("MST: could not open %s\n",filename);
		exit(1);
	}

   fscanf(distfile, "%d", &size);

	if (size<=0 || size>MAX_NUM_DIST) {
		printf ("size=%d in file %s not allowed\n",size,filename);
		exit(1);
	}

   for(i=0; i < size; i++) {
      for(j=0; j< size; j++) {
        fscanf(distfile, "%f", &temp);
		  /*
			w[i][j]= Pam_Estimated_From_Observed(temp);

			w[i][j]= temp;
		  */
			w[i][j]= temp*temp;
      }
   }

  fclose(distfile);

  return size;

}


/*
   This file is part of the Pratt program source code
   Copyright 1996 Inge Jonassen, Dept. of Informatics
   University of Bergen. 
   email: inge@ii.uib.no
   More information on Pratt: 
   http://www.ii.uib.no/~inge/Pratt.html
*/
