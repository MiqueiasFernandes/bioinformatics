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
#include <pthread.h>
#include <sys/sysinfo.h>
#include <stdbool.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <math.h>
#include "kseq.h"

#define USE_ALL_CORES false
#define MONTADOR 'L' //Illumina 1.8+ Phred+33
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

KSEQ_INIT(gzFile, gzread)

typedef struct seq {
	char nome[100];
	char seq[1000];
	char qal[1000];
	char primer1[100];
	char primer2[100];
	int tam;
	int f;
	int pos;
	char sequencia_nome[1000];
	char sequencia_seq[1000];
	char sequencia_qual[1000];
} SEQ;

typedef struct tdef {
	int n;
	long int inicio;
	long int fim;
	long int lidos;
	long int seqs;
	long int excedente;
	char* fastq;
	SEQ* sequences;
	int sequencesc;
	char* save;
	bool isAlive;
} TDEF;

int strContains(char* source, char* target) {
	int source_sz = strlen(source);
	int target_sz = strlen(target);
	if (source_sz < 1 || target_sz < 1) {
		return source_sz == target_sz ? 0 : -1;
	}
	bool match = false;
	int i, j;
	for (i = 0; i < (source_sz - (target_sz - 1)); ++i) {
		for (j = 0; j < target_sz; ++j) {
			if (!(match = (source[i + j] == target[j])))
				break;
		}
		if (match)
			break;
	}
	return match ? i : -1;
}

int includes(char* a, char* b) {
	return strlen(a) > strlen(b) ? strContains(a, b) : strContains(b, a);
}

void *thread(void *x_void_ptr) {

	TDEF* tdef = (TDEF *) x_void_ptr;
	int m;
	gzFile fp;
	kseq_t *fseq;
	long int l;
	fp = gzopen(tdef->fastq, "r");

	if (tdef->inicio != (l = gzseek(fp, tdef->inicio, SEEK_SET))) {
		printf("thread %d falhou ao pular para frente: %d.\n", tdef->n, l);
	} else {
		//printf("thread %d iniciar em %ld - esta em %ld.\n", tdef->n, tdef->inicio, gztell(fp));
		fseq = kseq_init(fp);
		while ((l = kseq_read(fseq)) >= 0) {
			tdef->seqs++;
			tdef->lidos = gztell(fp) - tdef->inicio;
			long int i, j, tmp;
			for (i = 0; i < tdef->sequencesc; ++i) {
				if ((tmp = strContains(fseq->seq.s, tdef->sequences[i].seq))
						>= 0) {
					if (tdef->sequences[i].f < 1)
						printf("sequencia %s encontrada.\n", tdef->sequences[i].nome);
					tdef->sequences[i].f++;
					for (j = 0; j < tdef->sequences[i].tam; ++j) {
						tdef->sequences[i].qal[j] = fseq->qual.s[tmp + j];
					}
					tdef->sequences[i].pos = tmp;
					tdef->sequences[i].qal[++j] = '\0';
					sprintf(tdef->sequences[i].sequencia_nome, "%s %s",
							fseq->name.s,
							fseq->comment.l ? fseq->comment.s : "");
					if (fseq->qual.l < 1000) {
						strcpy(tdef->sequences[i].sequencia_qual, fseq->qual.s);
						strcpy(tdef->sequences[i].sequencia_seq, fseq->seq.s);
					} else {
						for (m = 0; m < 900; m++) {
							tdef->sequences[i].sequencia_qual[m] = fseq->qual.s[m];
							tdef->sequences[i].sequencia_seq[m] = fseq->seq.s[m];
						}
						tdef->sequences[i].sequencia_qual[m] = '\0';
						tdef->sequences[i].sequencia_seq[m] = '\0';
						strcat(	tdef->sequences[i].sequencia_qual, " ... ");
						strcat(	tdef->sequences[i].sequencia_seq, " ... ");				
					}
					if (tdef->save) {
						FILE* f = fopen(tdef->save, "a");
						fprintf(f, ">%s\n%s\n%s\n%s\n%d VEZES ENC, POS %d\n", tdef->sequences[i].nome, tdef->sequences[i].sequencia_nome, fseq->seq.s, fseq->qual.s, tdef->sequences[i].f, tmp);						
						fclose(f);					
					}
				}
			}
			if (tdef->lidos >= (tdef->fim - tdef->inicio))
				break;
		}
		tdef->excedente = (tdef->inicio + tdef->lidos) - (tdef->fim + 1);
	}
	kseq_destroy(fseq);
	gzclose(fp);

	tdef->isAlive = false;
	return NULL;
}

long int getSizeFile(char* file) {
	FILE* fp = fopen(file, "r");
	fseek(fp, 0L, SEEK_END);
	long int res = 0;
	res = ftell(fp);
	fclose(fp);
	return res;
}

long int getSequences(char* file, long int size) {
	printf("Pre processando arquivo (%ld bytes)...\n%ld a processar...", size, size);
	gzFile fp;
	kseq_t *seq;
	fp = gzopen(file, "r");
	seq = kseq_init(fp);
	long int cont = 0;
	int l = 0;
	long int perc = 0;
	while ((l = kseq_read(seq)) >= 0){
		cont++;
		perc = gztell(fp);
		printf("\r%ld a processar [%ld] (%d%)...                ", size - perc, cont, (perc * 100) / size);
	}
	kseq_destroy(seq);
	gzclose(fp);
	printf("\n%ld seq encontradas.\n", cont);
	return cont;
}

long int ajustar(char* file, long int local) {

	FILE* fp = fopen(file, "r");
	fseek(fp, local, SEEK_SET);
	char buffer[100000];
	long int primeiralinha = 0;
	long int segundalinha = 0;
	int linhas = 0;
	long int skip = 0;
	long int ret = -1;
	bool arroba = false;

	while (fgets(buffer, 100000, fp) && (strcmp(buffer, "+\n") != 0)) { ///procura linha com +
		if (linhas == 0) {
			arroba = buffer[0] == '@';
			primeiralinha = strlen(buffer);
		} else if (linhas == 1) {
			segundalinha = strlen(buffer) + 1;
		}

		linhas++;
		skip += strlen(buffer);
	}

	if (strcmp(buffer, "+\n") == 0) { /// achou +
		skip += 2;

		switch (linhas) {
		case 0: /// pular apenas a proxima linha
		case 1:
			fgets(buffer, 100000, fp);
			ret = local + skip + strlen(buffer);
			break;
		case 2:
			if (arroba) {
				ret = local;
			} else {
				fgets(buffer, 100000, fp);
				ret = local + skip + strlen(buffer);
			}
			break;
		case 3:
			ret = local + primeiralinha;
			break;
		case 4:
			ret = local + segundalinha;
			break;
		}

	} else { /// nao achou +
		printf("falha no arquivo fastq\n");
	}

	fseek(fp, ret, SEEK_SET);
	fgets(buffer, 100000, fp);

	fclose(fp);
	return ret;
}

char* getSeq1(char* seq1) {
	int i, k =0;
	for(i = 0; i < strlen(seq1); i++)
		if (seq1[i] == 'A' || seq1[i] == 'C' || seq1[i] == 'T' || seq1[i] == 'G' || seq1[i] == 'N' || seq1[i] == '-') 
			seq1[k++] = seq1[i];
	seq1[k] = '\0';
	return seq1;
}

char* getPadding(char* b, int p) {
	b[p] = '\0';
	for (p -= 1; p >= 0; --p) {
		b[p] = ' ';
	}
	return b;
}

int procQ(char str) {
	int q = (int) str;
	switch(MONTADOR) {
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
		printf("Modo desconhecido \"%c\"\n", MONTADOR);
		return -1;
	}
	return q;
}

int processar(char* b1, char* b2, char* b3, SEQ* seq, int linha) {
	printf("\n%s%s%s", b1, b2, b3);
	int i, j, tmp, snips1 = 0, snips2 = 0;
	char s1[1000], s2[1000];
	char p[1000];
	char sp[1000];
	char bb1[1000];
	char bb2[1000];
	char f1[1000], t1[1000], q1[1000];
	char f2[1000], t2[1000], q2[1000];
	int p11[1000], p12[1000];
	int p21[1000], p22[1000];
	char pr1[100];
	char pr2[100];
	int strl1 = strlen(seq->primer1);
	int strl2 = strlen(seq->primer2);
	strcpy(pr1, seq->primer1);
	strcpy(pr2, seq->primer2);
	char seq1[1000];
	char seq2[1000];
	strcpy(seq1, b1);
	strcpy(seq2, b3);
	getSeq1(seq1);
	getSeq1(seq2);
	int pos = strContains(seq->sequencia_seq, seq->seq);
	getPadding(p, pos);
	int c, r = 0, m1 = -1, m2 = -1, q;
	char pd1[1000], pd2[1000];

	if ((tmp = strContains(b3, pr1)) >= 0) {
		j = strContains(seq->seq, pr1);
		getPadding(sp, j);
		sprintf(s1, "%s%s%s", p, sp, pr1);
		c = strlen(s1) - strl1;
		for (i = 0; i < strl1; ++i) {
			bb1[i] = seq->qal[j + i];
			if (b2[tmp + i] != '|') {
				f1[snips1] = b1[tmp + i];
				t1[snips1] = b3[tmp + i];
				q1[snips1] = seq->qal[j + i];
				p11[snips1] = c + i + 1;
				p12[snips1] = i + 1;
				snips1++;
				pr1[i] = s1[c + i] = 48 + snips1 % 10;
			}
		}
		bb1[i] = '\0';
		r++;
		m1 = strContains(seq2, seq->primer1);
		if (m1 >= 0) {
			getPadding(pd1, m1);		
		}
	} else sprintf(s1, "PRIMER 1 NÃO ENCONTRADO (%s)", pr1);

	if ((tmp = strContains(b3, pr2)) >= 0) {
		j = strContains(seq->seq, pr2);
		getPadding(sp, j);
		sprintf(s2, "%s%s%s", p, sp, pr2);
		c = strlen(s2) - strl2;
		for (i = 0; i < strl2; ++i) {
			bb2[i] = seq->qal[j + i];
			if (b2[tmp + i] != '|') {
				f2[snips2] = b1[tmp + i];
				t2[snips2] = b3[tmp + i];
				q2[snips2] = seq->qal[j + i];
				p21[snips2] = c + i + 1;
				p22[snips2] = i + 1;
				snips2++;
				pr2[i] = s2[c + i] = 48 + snips2 % 10;
			}
		}
		bb2[++i] = '\0';
		r++;
		m2 = strContains(seq2, seq->primer2);
		if (m2 >= 0) {
			getPadding(pd2, m2);		
		}
	} else sprintf(s2, "PRIMER 2 NÃO ENCONTRADO (%s)", pr2);
	
	printf(
			"\n>%s, encontrada %dX no arquivo fastq, L %d C %d\n%s\n%s\n%s%s\n%s\n%s\n%s\n",
			seq->nome, seq->f, linha, pos + 1, seq->sequencia_nome,
			seq->sequencia_seq, p, seq->seq, s1, s2, seq->sequencia_qual);
	if (m1 >= 0) {
		printf("   ***   primer 1   ***\n%s\n%s%s\n%s\nSEQ: %s\n", seq1, pd1, seq->primer1, seq2, seq->primer1);
		j = 0;
		printf("SNP: ");
		for (i = 0; i < strlen(pr1); ++i) {
			if (pr1[i] < 48 || pr1[i] > 58) {
				printf(" ");
			} else {
				printf("%c", f1[j++]);
			}
		}
		printf("\nQUA: %s\n", bb1);
		for (i = 0; i < snips1; ++i) {
			q = procQ(q1[i]);
			printf("SNP %d : %c -> %c = %c (%d ~%d%) [offset seq: %d prim: %d]\n", i + 1,
					f1[i], t1[i], q1[i], q, ((q * 100) / 41), p11[i], p12[i]);
		}
	}
	if (m2 >= 0) {
		printf("   ***   primer 2  %d ***\n%s\n%s%s\n%s\nSEQ: %s\n", m2, seq1, pd2, seq->primer2, seq2, seq->primer2);
		j = 0;
		printf("SNP: ");
		for (i = 0; i < strlen(pr2); ++i) {
			if (pr2[i] < 48 || pr2[i] > 58) {
				printf(" ");
			} else {
				printf("%c", f2[j++]);
			}
		}
		printf("\nQUA: %s\n", bb2);
		for (i = 0; i < snips2; ++i) {
			q = procQ(q2[i]);
			printf("SNP %d : %c -> %c = %c (%d ~%d) [offset seq: %d prim: %d]\n", i + 1,
					f2[i], t2[i], q2[i], q, ((q * 100) / 41), p21[i], p22[i]);
		}
	}
	printf("     ####          ####          ####          ####          ####          ####     \n\n");
	return r > 0 ? 1 : 0;
}

//https://stackoverflow.com/questions/29787310/does-pow-work-for-int-data-type-in-c
long int int_pow(int base, int exp)
{
    long int result = 1;
    while (exp)
    {
        if (exp & 1)
           result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

long int getInt(char* str) {
	int i, j;
	long int num = 0;
	for(i = strlen(str) - 1, j = 0; i >= 0; --i, j++) {
		num += (int_pow(10, j) * (str[i] - '0')); 
	}
	return num;
}

int principal(int argc, char *argv[], int vez) {
	int r = 0;
	printf("------Buscar qualidade da sequencia S (%dX)------\n", vez);
	if (argc > 3) {

		//PADRAO PARA ENTRADA:
		//L1) #gene_nome
		//L2) ACTGCG ------> seq maior
		//L3) ACT  --------> primer1
		//L4) ACT  --------> primer2
		if (argc > 6) printf("sequencias limitadas a %ld.\n", getInt(argv[6]));
		FILE* fseq = fopen(argv[2], "r");
		SEQ seqs[1000];
		int cont = 0;
		int max = 0;
		while (fscanf(fseq, "#%s\n%s\n%s\n%s\n", seqs[cont].nome,
				seqs[cont].seq, seqs[cont].primer1, seqs[cont].primer2) > 0) {
			if (argc > 6) {
				seqs[cont].seq[max = getInt(argv[6])] = '\0';
			} else {
				max = MAX(max, strlen(seqs[cont].seq));		
			}
			printf("%d encontrar: %s\n", cont + 1, seqs[cont].nome);
			seqs[cont].f = 0;
			seqs[cont].tam = strlen(seqs[cont].seq);
			cont++;
		}
		fclose(fseq);

		int cr = get_nprocs();
		int ts = cr - (USE_ALL_CORES ? 0 : 1);
		long int fileSize = getSizeFile(argv[1]);
		long int sizePorThread = fileSize;
		long int nfastseqs = (argc > 4) ? getInt(argv[4]) : getSequences(argv[1], fileSize);

		if (nfastseqs > 10 && ts > 1) {
			sizePorThread = fileSize / ts;
		}

		printf("sequencias no fastq: %ld\nsequencias a encontrar: %d\nnucleos a usar:  %d de %d\ntamanho do arquivo: %ld frag: %ld\n",nfastseqs, cont, ts, cr, fileSize, sizePorThread);

		pthread_t threads[ts];
		TDEF* tdefs[ts];
		int i;
		long int tmp = 0;

		printf("iniciando threads: 1");
		for (i = 0; i < ts; ++i) {

			TDEF *tdef = tdefs[i] = malloc(sizeof(TDEF));

			tdef->n = i + 1;
			tdef->inicio = ajustar(argv[1], tmp);
			tdef->fim = tmp + sizePorThread - 1;
			tdef->lidos = 0;
			tdef->seqs = 0;
			tdef->fastq = argv[1];
			tdef->sequences = seqs;
			tdef->sequencesc = cont;
			tdef->save = (argc > 5) ? argv[5] : NULL;
			tdef->excedente = 0;
			tdef->isAlive = true;
			if (pthread_create(&threads[i], NULL, thread, tdef)) {
				printf("Error creating thread %d\n", tdef->n);
				return 1;
			} else {
				if (i > 0) printf(", %d", tdef->n);
				tmp += sizePorThread;
			}
		}
		printf("\n%d threads iniciadas\n", i);

		long int lido = 0;
		int term = 0;
		//long int exc = 0;
		long int proc = 0;
		int per = -1;
		int enc = 0;
		while (lido < fileSize && term < ts) {
			lido = proc = term = enc = 0;
			for (i = 0; i < ts; ++i) {
				lido += tdefs[i]->lidos;
				proc += tdefs[i]->seqs;
				///exc += tdefs[i]->excedente;
				if (!tdefs[i]->isAlive)
					term++;
			}
			for (i = 0; i < cont; i++) {
				enc += seqs[i].f > 0 ? 1 : 0;			
			}
			if (((lido * 100) / fileSize) > per) 
				printf("PER %d%c THR: %d. SEQ: %ld. ENC: %d\n", (per = ((lido * 100) / fileSize)),'%', ts -term, nfastseqs -proc, enc);
			sleep(1);
		}
		printf("aguardando threads terminarem...\n");

		for (i = 0; i < ts; ++i) {
			pthread_join(threads[i], NULL);
		}

		for (i = 0; i < cont; i++) {
			if (seqs[i].f < 1) {
				printf("SEQ: %s : NÃO FOI ENCONTRADA\n", seqs[i].nome);	
			} else {
				printf("SEQ: %s : encontrada %dx\n", seqs[i].nome, seqs[i].f);		
			}		
		}
		printf("\nbuscando sequencias no arquivo %s...\n", argv[3]);
		//encontar arquivo qual espec
		FILE* fal = fopen(argv[3], "r");
		char b1[1000];
		char b2[1000];
		char b3[1000];
		long int k, linha = 3;
		bool eh = true;
		if (fgets(b1, 1000, fal)) {
			if (fgets(b2, 1000, fal)) {

				while (fgets(b3, 1000, fal)) {
				eh = true;
				for (i = 0; i < strlen(b2); ++i)
					if (!(b2[i] == '|' || b2[i] == ' ' || b2[i] == '\n')) 
						eh = false;

					for (i = 0; eh && i < cont; ++i) {
						printf("\rprocessando seq %s para linha %d.                  ", seqs[i].nome, linha);
						if (seqs[i].f > 0) {
							if (MAX(strContains(b3, seqs[i].primer1), strContains(b3, seqs[i].primer1)) >= 0) {
								r += processar(b1, b2, b3, &seqs[i], linha);
							}
						}
					}

					strcpy(b1, b2);
					strcpy(b2, b3);
					linha++;
				}
			}
		}
		fclose(fal);
		printf("\n%d sequencias encontradas com sucesso na etapa 1 e 2.\n", r);
		r = cont - r;
	} else {
		printf(
				"USAGE: ./qualityfinder  [path fastq] [path sequences] [path alinhado] [path save rep] [num limit seq] [num iteração] [num decremento iterativo].\n");
	}
	return r;
}

int main(int argc, char *argv[]) {
	if (argc > 7) {
		int it = argc > 8 ? getInt(argv[8]) : 1;
		printf("modo iterativo limitado a %ldX (%d).\n", getInt(argv[7]), it);
		int i;		
		for (i = 0; i < getInt(argv[7]); i++) {
			if (principal(argc, argv, i + 1) < 1){
				break;			
			}
			sprintf(argv[6], "%d", getInt(argv[6]) - it);
		}
	} else {
		principal(argc, argv, 1);
	}	
	printf("terminado com sucesso!\nby mikeias.net\n");
	return 0;

}

