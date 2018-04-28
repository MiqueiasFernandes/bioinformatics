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
 Name        : qualityfinder.c
 Author      : Miquéias Fernandes
 Version     : 1.0 25/04/2018
 Description : Get quality of sequences in .fastq
 Compile     : gcc -g -O2 -lz -pthread qualityfinder.c
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <stdbool.h>

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef struct p {
	char id[100];
	char seq1[100];
	char seq2[100];
	int enc;
} PRIMER;

int strContains(char* source, char* target) {
	int source_sz = strlen(source);
	int target_sz = strlen(target);
	if (source_sz < 1 || target_sz < 1) {
		return source_sz == target_sz ? 0 : -1;
	}
	bool match = false;
	int i, j;
	for (i = 0; i < (source_sz - (target_sz - 1)); ++i) {
		for (j = 0; j < target_sz; ++j)
			if (!(match = (source[i + j] == target[j])))
				break;
		if (match)
			break;
	}
	return match ? i : -1;
}

long int getSizeFile(char* file) {
	FILE* fp = fopen(file, "r");
	fseek(fp, 0L, SEEK_END);
	long int res = 0;
	res = ftell(fp);
	fclose(fp);
	return res;
}

int main(int argc, char *argv[]) {

	printf("------Buscar GENE em arquivo FASTA ------\n");
	int k = 0, t1, t2, sqc = 0, ac = 0, i;
	char id[1000];
	char seq[1000000];
	char str[1000];

	char b1[100];
	char b2[100];
	char b3[100];
	char b4[100];

	PRIMER ps[100];
	FILE* f1 = fopen(argv[1], "r");
	while (fscanf(f1, ">%s %s %s\n%s\n>%s %s %s\n%s\n>%s %s %s %s\n%s\n", b4,
			b1, b4, b2, b4, b4, b4, b4, b4, b4, b4, b4, b3) > 0) {
		strcpy(ps[sqc].id, b1);
		strcpy(ps[sqc].seq1, b2);
		strcpy(ps[sqc].seq2, b3);
		ps[sqc].enc = 0;
		sqc++;
	}

	fclose(f1);

	printf("Foram encontrados [%d] seq no arquivo\n", sqc);

	for (k = 0; k < sqc; k++) {
		printf("%s, ", ps[k].id);
	}
	printf("\n");

	int a;
	FILE* fseq;
	long int cont, size;
	int de, ate;
	FILE* saida = fopen(argv[argc - 1], "w");
	char* r;
	for (a = 2; a < (argc - 1); a++) {
		printf("\nprocessando arquivo %s\n", argv[a]);
		size = getSizeFile(argv[a]);
		fseq = fopen(argv[a], "r");
		cont = 0;
		str[0] = seq[0] = '\0';
		while ((r = fgets(str, 1000, fseq)) > 0) {
			cont += strlen(str);
			if (str[0] == '>' || cont == size) {
				if (str[0] != '>') {
					if (str[strlen(str) - 1] == '\n') str[strlen(str) - 1] = '\0';
					strcat(seq, str);
				}
				if (strlen(id) > 0 && strlen(seq) > 0) {
					for (k = 0; k < sqc; k++) {
						if (((t1 = strContains(seq, ps[k].seq1)) >= 0) && ((t2 =
								strContains(seq, ps[k].seq2)) >= 0)) {
							printf("\n[%dX] %s => %s\n", ++ps[k].enc, ps[k].id, id);
							if (ps[k].enc == 1)
								ac++;
							de = MIN(t1, t2);
							ate = MAX(t1, t2);
							printf("#SEQ ENTRE %d e %d (%dpb) [%d]\n", de, ate,
									ate - de, strlen(seq));
							fprintf(saida, ">%s %s %d e %d (%dpb)\n", ps[k].id,
									argv[a], de, ate, ate - de);
							for (i = de; i < ate; i++) {
								//printf("%c", seq[i]);
								fputc(seq[i], saida);
							}
							//printf("%s\n", ate == t1 ? ps[k].seq1 : ps[k].seq2);
							fprintf(saida, "%s\n",
									ate == t1 ? ps[k].seq1 : ps[k].seq2);
						}
					}
				}
				strcpy(id, str);
				id[strlen(id) - 1] = '\0';
				seq[0] = '\0';
			} else {
				str[strlen(str) - 1] = '\0';
				strcat(seq, str);
			}

			printf("\rprocessando %ld%c", ((cont * 100) / size), '%');
		}

		fclose(fseq);

	}
	fclose(saida);
	printf("\n\nRESUMO:\nENCONTRADOS: %d\n", ac);

	for (k = 0; k < sqc; k++) {
		printf("%s => %dX\n", ps[k].id, ps[k].enc);
	}

	printf("terminado com sucesso!\nby mikeias.net\n");
	return 0;

}

