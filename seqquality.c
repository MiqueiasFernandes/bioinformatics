/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 ============================================================================
 Name        : seqquality.c
 Author      : Miquéias Fernandes
 Version     : 1.0 25/04/2018
 Description : Get quality of sequences
 ============================================================================
*/

#include <stdio.h>
#include <string.h>

// https://en.wikipedia.org/wiki/FASTQ_format 
//  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
//  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
//  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
//  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ.....................
//  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
//  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
//  |                         |    |        |                              |                     |
// 33                        59   64       73                            104                   126
//  0........................26...31.......40                                
//                           -5....0........9.............................40 
//                                 0........9.............................40 
//                                    3.....9..............................41 
//  0.2......................26...31........41                              

#define HELP "\
 S - Sanger        Phred+33,  ! a I, raw reads typically (0, 40)\n\
 X - Solexa        Solexa+64, ; a h, raw reads typically (-5, 40)\n\
 I - Illumina 1.3+ Phred+64,  @ a h, raw reads typically (0, 40)\n\
 J - Illumina 1.5+ Phred+64,  C a h, raw reads typically (3, 41)\n\
 L - Illumina 1.8+ Phred+33,  ! a I, raw reads typically (0, 41)\n"

int processar(char* str, char modo, int ver) {
	if (ver > 0) printf("PROCESSANDO STRING: %s\n", str);
	int i, q;
	for (i=0; i < strlen(str); i++) {
		q = (int) str[i];
		switch(modo) {
			case 'S': 
			case 's': 
			case 'L': 
			case 'l':
			q -= (int) '!';
			break;
			case 'X': 
			case 'x':
			case 'I': 
			case 'i':  
			case 'J': 
			case 'j': 
			q -= (int) '@'; 
			break;
			default:
			printf("Modo desconhecido \"%c\"\n", modo);
			return -1;
		}
		if (ver > 0) printf("[%i] : %c => %d [~%d%]\n", i, str[i], q, (q * 100) / 41);
	}
	return 0;
}

int main(int argc, char* argv[]) {
	if (argc < 2 || processar("test", argv[1][0], 0) < 0) 
		return printf("usage ./qual [tipo] [seq]*\nOnde tipo pode ser:\n%s", HELP);
	int i;
	if (argc > 2)
		for ( i = 2; i < argc; i++)
			processar(argv[i], argv[1][0], 1);
	char str[1000];
	do{
		fflush(stdout);
		fflush(stdin);
		printf("Digite uma string ou q para sair...\n");
		scanf("%s", str);
		if (str[0] == 'q') return 0;
		processar(str, argv[1][0], 1);
	} while(1);
}
