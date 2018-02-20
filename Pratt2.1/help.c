#include <stdio.h>
#include "help.h"

typedef struct 
{
	char *subject;
	char *content;
} t_Help;

#define NUMBER_HELP 12
t_Help Help_Var[NUMBER_HELP]= 
{ 
{"HELP", "\nGeneral help on Pratt\n"
"=====================\n\n"
"Pratt is a tool for discovering patterns that match some minimum number\n"
"of a set of sequences. The sequences are input all in one file in either\n"
"FastA or Swiss-Prot format.\n\n"
"You specify:\n"
"  The minimum number of sequences to match a pattern (using options M or %)\n"
"  Restrictions on patterns to be found (using options BN,PL,PC,PX,FN,FL,FP options)\n"
"  Greedyness of the search (using the E option).\n\n"
"The patterns found by Pratt will be output to a file.\n"
"You can set the name of this file using option OF.\n\n"

"You can also set\n"
"  the number of patterns to be reported (option ON)\n"
"  the format of the patterns (OP)\n"
"  the number of patterns for which an alignment is shown (OA).\n\n"

"Pratt can also output a summary of where in the sequences each of\n"
"the patterns match. This will be done when option M is switched on.\n"
"Options MV and MR allows you to adjust the output.\n\n"

"Optionally, Pratt can be instructed only to look for patterns that\n"
"  are consistent with an alignment, or match a special \n"
"  query sequence (option J). \n\n"

"References: \n\n"

"  I.Jonassen, J.F.Collins, D.G.Higgins.\n"
"  Protein Science 1995;4(8):1587-1595.\n\n"

"  I. Jonassen \n"
"  Efficient discovery of conserved patterns using a pattern graph\n"
"  Submitted to CABIOS\n\n"

"For more information on the use of any of the options in the menu,\n"
"type 'help' <option>.\n"
},
{"C", 
"The C parameters control how many sequences a pattern should match to\n"
"be considered by Pratt.\n"

"\nCM: sets the minimum number of sequences to match a pattern.\n\n"
"Pratt will only report patterns that match at least the chosen\n"
"number of the sequences that you have input. Pratt will not allow\n"
"you to choose a value higher than the number of sequences input.\n\n"

"\nC%: sets the minimum percentage of the input sequences that should match\n"
"a pattern.\n\n"
"If you set this to, say 80, Pratt will only report patterns \n"
"matching at least 80 % of the sequences input. \n\n"
},
{"G",
"Allows the use of an alignment or a query sequence to restrict the pattern search.\n"
"\n"
"If G is set to al or query, another option GF will appear\n"
"allowing the user to give the name of a file containing a\n"
"multiple sequence alignment (in Clustal W format), or a \n"
"query sequence in FastA format (without annotation).n"
"\n"
"Only patterns consistent with the alignment/matching the query\n"
"sequence will be considered.\n"
"\n"
"Loosely a pattern is considered consistent with the alignment if \n"
"- each symbol in the pattern (e.g. A) corresponds to a ungapped column\n"
"  in the alignment where all the characters match the pattern symbol\n"
"  (in the example, A).\n"
"- the wildcards in the pattern are compatible with the number of\n"
"  residues between the corresponding columns in the alignment.\n"
"\n"
"For instance the pattern A-x(2,3)-B is consistent with the alignment\n"
"   ALVGB\n"
"   AG-LB\n"
"   ALD-B\n"
"\n"
"For more details see\n"
"\n"
"  I. Jonassen \n"
"  Efficient discovery of conserved patterns using a pattern graph\n"
"  Submitted to CABIOS\n"
},
{"B",
"\nUsing the B options (BN,BI,BF) on the menu you can control which \n"
"pattern symbols will be used during the initial pattern\n"
"search and during the refinement phase. \n"
"In the pattern C-x(2)-[DE], C and [DE] are the symbols.\n"
"\n"
"The pattern symbols that can be used, are read from a file \n"
"if the BI option is set, otherwise a default set will be used.\n"
"\n"
"The default set has as the 20 first elements, the single amino acid\n"
"symbols, and it also contains a set of ambiguous symbols, each \n"
"containing amino acids that share some physio-chemical properties\n"
"\n"
"If BI is set, option BF will appear to allow you to give the\n"
"name of the file.  In the file each symbol is given on a separate\n"
"line contatining the letters that the symbol should match. \n"
"For instance the file could be:\n"
"\n"
"C\n"
"DE\n"
"\n"
"and only patterns with the symbols C and [DE] would\n"
"be considered. During the initial search, pattern symbols\n"
"corresponding to the first BN lines can be used.\n"
"Increasing BN will slow down the search and increase the\n"
"memory usage, but allow more ambiguous pattern symbols.\n"
},
{"P",
"The P options are for controlling the patterns to be considered\n"
"by Pratt. See also the F options for controlling flexibility\n"

"\nOption PL: allows you to set the maximum length of a pattern. \n"
"The length of the pattern C-x(2,4)-[DE] is 1+4+1=6.\n"
"The memory requirement of Pratt depends on L; a higher L\n"
"value gives higher memory requirement.\n"

"\nOption PN: using this you can set the maximum number of symbols\n"
"in a pattern. The pattern C-x(2,4)-[DE] has 2 symbols (C and [DE]).\n"
"When PN is increased, Pratt will require more memory.\n"

"\nOption PX: Using this option you can set the maximum length of a \n"
"wildcard. Examples of wildcards and lengths are\n"
"   x       -  1\n"
"  x(10)    - 10\n"
" x(3,4)    -  4\n"
"\n"
"Increasing PX will increase the time used by Pratt,\n"
"and also slightly the memory required.\n"
},
{"F",
"The F options control flexible wildcards in the patterns.\n"

"\nOption FN: Using this option you can set the maximum number of\n"
"flexible wildcards (matching a variable number of\n"
"arbitrary sequence symbols). For instance x(2,4)\n"
"is a flexible wildcard, and the pattern\n"
"\n"
"   C-x(2,4)-[DE]-x(10)-F \n"
"\n"
"contains one flexible wildcard.\n"
"\n"
"Increasing FN will increase the time used by Pratt.\n"
"\nOption FL: you can set the maximum flexibility of\n"
"a flexible wildcard (matching a variable number of\n"
"arbitrary sequence symbols). For instance x(2,4) and\n"
"x(10,12) has flexibility 2, and x(10) has flexibility 0.\n"
"\n"
"Increasing FL will increase the time used by Pratt.\n"

"\nOption FP: Using option P you can set an upper limit on the product\n"
"of a flexibilities for a pattern. This is related to \n"
"the memory requirements of the search, and increasing\n"
"the limit, increases the memory usage.\n"
"\n"
"Some patterns and the corresponding product of flexibilities:\n"
"\n"
"C-x(2,4)-[DE]-x(10)-F      - (4-2+1)*(10-10+1)= 3\n"
"C-x(2,4)-[DE]-x(10-14)-F   - (4-2+1)*(14-10+1)= 3*5= 15\n"
},
{"E",
"Using the E parameter you can adjust the greediness of the\n"
"search. Setting E to 0 (zero), the search will be exhaustive.\n"
"Increasing E increases the greediness, and decreases the time\n"
"used in the search.\n"
},
{"R",
"When the R option is switched on, patterns found during the\n"
"initial pattern search are input to a refinement algorithm\n"
"where more ambiguous pattern symbols can be added.\n"
"\n"
"For instance the pattern \n"
"    C-x(4)-D \n"
"might be refined to\n"
"C-x-[ILV]-x-D-x(3)-[DEF]\n"
"\n"
"If the RG option is switched on, then ambiguous symbols listed in the\n"
"symbols file (or in the default symbol set -- see help for option B),\n"
"are used. If RG is off, only the letters needed to match the input\n"
"sequences are inlcuded in the ambiguous pattern positions.\n"
"\n"
"For example, if [ILV] is a listed allowed symbol, and [IL] is not, [IL] can\n"
"be included in a pattern if RG is off, but if RG is on, the full symbol\n"
"[ILV] will be included instead.\n"
},
{"O",
"The O options allow you to control the output from Pratt.\n"

"\nOption OF: allows you to specify the name of the file to \n"
"which Pratt will write its output\n"

"\nOption OP: when switched on, patterns will be output in \n"
"PROSITE style (for instance C-x(2,4)-[DE]). When switched\n"
"off, patterns are output in a simpler consensus pattern\n"
"style (for instance Cxx--[DE] where x matches exactly one\n"
"arbitrary sequence symbol and - matches zero or one arbitrary\n"
"sequence symbol).\n"

"\nOption ON: sets the max. nr of patterns to be found by Pratt.\n"

"\nOption OA: sets the max. nr of patterns for which Pratt is to\n"
"produce an alignment of the sequence segments matching it.\n"
},
{"M",
"\nIf the M option is set, then Pratt will print out the location \n"
"of the sequence segments matching each of the (maximum 52) best \n"
"patterns. The patterns are given labels A, B,...Z,a,b,...z in order\n"
"of decreasing pattern score. Each sequence is printed on a line, \n"
"one character per K-tuple in the sequence. If pattern with label C \n"
"matches the third K-tuple in a sequence C is printed out. If several \n"
"patterns match in the same K-tuple, only the best will be printed. \n"

"\nOption MR: sets the K value (ratio) used for printing the summary\n"
"information about where in each sequence the pattern matches are found.\n"

"\nOption MV: if set, the output is printed vertically instead of horizontally,\n"
"vertical output can be better for large sequence sets.\n"
},
{"S",
"\nThe S option allows you to control the scoring of patterns.\n"
"\nThere are five possible scoring schemes to be used:\n"
"info - patterns are scored by their information content as defined in\n"
"       (Jonassen et al, 1995). Note that a pattern's score is independent\n"
"       of which sequences it matches.\n"
"mdl  - patterns are scored by a Minimum Description Length principle\n"
"       derived scoring scheme, which is related to the one above, but\n"
"       penalises patterns scoring few sequences vs. patterns scoring many.\n"
"       Parameters Z0-Z3 appears when this scoring scheme is used.\n"
"tree - a pattern is scored higher if it contains more information and/or\n"
"       if it matches more diverse sequences. The sequence diversity is\n"
"       calculated from a dendrogram which has to be input.\n"
"dist - similar to the tree scoring, except a matrix with pairwise the\n"
"       similarity between all pairs of input sequences are used instead\n"
"       of the tree. The matrix has to be input."
"ppv  - a measure of Positive Predictive Value - it is assumed that the\n"
"       input sequences consitute a family, and are all contained in the\n"
"       Swiss-Prot database. PPV measures how certain one can be that a\n"
"       sequence belongs to the family given that it matches the pattern.\n"
"For the last three scoring schemes, an input file is needed and option SF\n"
"appears allowing the user to set his own file name.\n"
},
{"OPTIONS",
"\n All parameters that can be set using Pratt's menu can also be set\n"
"using options on the commandline. The general format is \n"
"\n"
"    -<menu description> <value>\n"
"\n"
"where the \n"
"\n"
" - menu description is the letter(s) used to set\n"
"   a parameter when using the menu, and\n"
"\n"
" - value is the value you want to assign to the parameter.\n"
"\n"
"For on/off parameters, <value> should be either 'on' or 'off',\n"
"for filename parameters, <value> should be a filename, and for\n"
"integer/real parameters, <value> should be a value.\n"
"\n"
"By giving the argument '-menu' you will also be presented with\n"
"the menu. Otherwise when command line options are given, the\n"
"menu will not be shown.\n"
}};



void general_help()
{
	int i;
	for (i=0;i<NUMBER_HELP;i++) {
		if (strcmp("OPTIONS",Help_Var[i].subject)==0) {
			printf ("%s\n",Help_Var[i].content);
		}
	}
}

void specific_help(char *subject)
{
	int i;
	char buf[50];

	for (i=0;i<NUMBER_HELP;i++) {
	/*
		if (strcmp(subject,Help_Var[i].subject)==0) {
	*/
		if (subject[0]==Help_Var[i].subject[0]) {
			printf ("Help on %s: \n\n%s\n\n",subject,Help_Var[i].content);

			printf ("Hit return to return to the menu.\n");

			fgets(buf,50,stdin);
			return;
		}
	}
	printf ("\nSorry, no help available on %s\n\n",subject);
	printf ("Hit return to return to the menu.\n");

	fgets(buf,50,stdin);
}
