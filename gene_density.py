
class Gene:
    def __init__(self, line):
        self.raw = line.split('\t')
        self.chr = self.raw[0]
        self.start = int(self.raw[1])
        self.end = int(self.raw[2])
    def contem(self, pos):
        return self.start <= pos and self.end >= pos

chrs = {}

for l in open('../genes_pos'):
    if l.startswith('Chr'):
            gene = Gene(l)
            if not gene.chr in chrs:
                chrs[gene.chr] = [gene]
            else:
                chrs[gene.chr].append(gene)

file = open("gene_density", "w")

for chrm in chrs:
    genes = chrs[chrm]
    genes.sort(key=lambda gene: gene.start)
    ini = min([g.start for g in genes])
    end = max([g.end for g in genes])
    print("%s %d => %d - %d" % (chrm, len(genes), ini, end), end=" ")
    numeros = {}
    for gene in genes:
        
        maior_q_ini = list(filter(lambda x: x >= gene.start, numeros))
        maior_q_end = list(filter(lambda x: x >= gene.end, numeros))
        
        if len(maior_q_ini) < 1:
            numeros[gene.start] = 1
            numeros[gene.end] = 1
            continue
        
        if len(maior_q_end) < 1:
            numeros[gene.end] = 1
            continue
        maior_q_ini.extend(maior_q_end)
        for num in set(maior_q_ini):
            numeros[num] += 1
    
    num_lim = {}
    num = 0
    maxi = 0
    for n in numeros:
        if numeros[n] != num:
            num = numeros[n]
            num_lim[n] = num
            maxi = max(maxi , num)
    
    for reg in num_lim:
        file.write("%s\t%d\t%d\t%.2f\n" % (chrm, reg, reg, num_lim[reg] / maxi))
    print(" (max %d overlap)" % maxi)
                
file.close()       


