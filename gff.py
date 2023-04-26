import sys

class GFF:
  def __init__(self, file):
    self.file = file
    self.valid = None
    self.invalid = []
    self.extras = []
    self.dup = []
    self.gene2mrna = []
    self.ptnas = None
    self.mrna2ptn = {}
    self.mrna2gene = {}
    self.genes = 0
    self.cgenes = 0
    self.vgenes = 0

  def valid_lines(self,):
    if not self.ptnas is None:
      with open(self.ptnas) as r:
        self.ptnas = set([l.strip() for l in r.readlines() if len(l) > 4])
        print('Qtd Ptnas:', len(self.ptnas))
    ## https://www.ensembl.org/info/website/upload/gff3.html
    cont = 0
    with open(self.file) as f:
      for l in f.readlines():
        cont += 1

        if l.startswith('#'):
          continue
        
        if l.count('\t') != 8:
          yield False, cont, 'few then 9 fields cont'
          continue

        fields = l.strip().split('\t')

        if (not fields[2] == 'gene') and (
            not fields[2] == 'mRNA') and (
            not fields[2] == 'exon') and (
            not fields[2] == 'CDS'):
            self.extras.append(l)
            continue

        ## field 1 must have a chr name
        if len(fields[0].strip()) < 1:
          yield False, cont, 'field CHR error ' + fields[0]
          continue

        ## field 4, 5 must have positive integer
        try:
          start, end = int(fields[3]), int(fields[4])
          assert start > 0 and end > 0
        except:
          yield False, cont, f'field not positive integer: {fields[3]} {fields[4]}'
          continue

        ## field 7 must have + or -
        strand = fields[6].strip()
        if (not strand == '+') and (not strand == '-'):
          yield False, cont, f'field 7 (strand) error: {fields[6]}, valid: +,-'
          continue

        ## field 9 must have ID 
        if (not fields[8].startswith('ID=')) and (not ';ID=' in fields[8]):
          yield False, cont, 'error not ID in field 9'
          continue

        ## field 9 must have Parent if is mRNA, exon or CDS 
        if fields[2] == 'mRNA' or fields[2] == 'exon' or fields[2] == 'CDS':
          if (not ';Parent=' in fields[8]) and (not fields[8].startswith('Parent=')):
            yield False, cont, 'field 9 without Parent:' + fields[8]
            continue
        yield True, cont, fields

  def valid_genes(self):
    know_genes = []
    know_mrnas = []
    know_mrnas_exon = []
    know_mrnas_cds = []
    invalid_mrnas = []
    mrna2gene = {}
    self.lines = [x for x in self.valid_lines()]
    ids, duplicated = [], []
    for s, k, line in self.lines:
      if not s:
        self.invalid.append(f'[{k}] => {line}\n')
        continue
      feature, anot = line[2], dict([x.split('=') for x in line[8].split(';')])
      if feature == 'gene' or feature == 'mRNA':
        if anot['ID'] in ids:
          duplicated.append(anot['ID'])
          self.dup.append(anot['ID'])
          continue
        else:
          ids.append(anot['ID'])
      if feature == 'gene':
        know_genes.append(anot['ID'])
        self.genes += 1
      if feature == 'mRNA':
        know_mrnas.append(anot['ID'])
        mrna2gene[anot['ID']] = anot['Parent']
      if feature == 'exon':
        know_mrnas_exon.append(anot['Parent'])
      if feature == 'CDS':
        know_mrnas_cds.append(anot['Parent'])
      if not self.ptnas is None and 'protein_id' in anot:
        if not anot['protein_id'] in self.ptnas:
          invalid_mrnas.append(anot['Parent'])
        else:
          self.mrna2ptn[anot['Parent']] = anot['protein_id']
    
    if len(ids) != len(set(ids)):
      raise Exception('Multiple feature to one ID.')
    
    ## valid if mRNA has at 1 exon & 1 CDS (coding mRNA)
    valid_mrnas = set(know_mrnas_exon).intersection(set(know_mrnas_cds))

    ## valid if mRNA has own line
    valid_mrnas.intersection_update(set(know_mrnas))
    
    ## valid if gene has own line
    ms = []
    for m, g in mrna2gene.items():
      if g in know_genes:
        ms.append(m)
    valid_mrnas.intersection_update(set(ms))

    ## valid if not ID duplicated
    valid_mrnas.difference_update(set(duplicated))

    ## valid if has ptnas
    valid_mrnas.difference_update(set(invalid_mrnas))

    ## valid if gene has at 1 valid mRNA
    valid_genes = set([mrna2gene[x] for x in valid_mrnas])

    dup_genes = valid_genes.intersection(set(duplicated))
    mrnaerr = [m for m, g in mrna2gene.items() if g in dup_genes]
    valid_genes.difference_update(dup_genes)
    valid_mrnas.difference_update(set(mrnaerr))
  
    self.cgenes = len(set([mrna2gene[x] for x in set(know_mrnas_cds) if x in mrna2gene]))
    self.vgenes = len(valid_genes)

    for s, _, line in self.lines:
      if not s:
        continue
      feature, anot = line[2], dict([x.split('=') for x in line[8].split(';')])
      if feature == 'gene':
        if anot['ID'] in valid_genes:
          yield line
      elif feature == 'mRNA':
        if anot['ID'] in valid_mrnas:
          self.gene2mrna.append((anot['Parent'], anot['ID']))
          self.mrna2gene[anot['ID']] = anot['Parent']
          yield line
      elif anot['Parent'] in valid_mrnas: ## exon or CDS only here
        yield line

  def load(self):
    if self.valid is None:
      self.valid = [x for x in self.valid_genes()]
    return not self.valid is None and len(self.valid) > 2
  
  def __try_load(self):
    if not self.load():
      raise Exception('Unable to load gff ' + self.file)

  def stats(self):
    self.__try_load()
    print('All Genes', self.genes)
    print('Coding Genes', self.cgenes)
    print('Genes', self.vgenes)
    print('mRNAs', len([x for x in self.valid if x[2] == 'mRNA']))
    print('Exons', len([x for x in self.valid if x[2] == 'exon']))
    print('CDSs', len([x for x in self.valid if x[2] == 'CDS']))
    print('Extra lines', len(self.extras))
    print('Invalid lines', len(self.invalid))
    print('Duplicated IDS', len(self.dup))
  
  def persist(self, name):
    self.__try_load()
    self.name = name
    with open(f'{name}_clean.gff', 'w') as fw:
      fw.writelines(['\t'.join(l) + '\n' for l in self.valid])
    with open(f'{name}_invalid.txt', 'w') as fw:
      fw.writelines(self.invalid)
    with open(f'{name}_extras.gff', 'w') as fw:
      fw.writelines(self.extras)
    with open(f'{name}_multiple.txt', 'w') as fw:
      fw.writelines([x+'\n' for x in self.dup])
    with open(f'{name}_gene2mrna.txt', 'w') as fw:
      fw.writelines([f'{x[0]}\t{x[1]}\n' for x in self.gene2mrna])
    if len(self.mrna2ptn) > 0:
      with open(f'{name}_gene2mrna2ptna.txt', 'w') as fw:
        fw.writelines([f'{x[0]}\t{x[1]}\t{self.mrna2ptn[x[1]]}\n' for x in self.gene2mrna])

  def __strand(self, seq, strand):
    if strand == '+':
      return seq
    assert seq.isalpha()
    a = {'A':'1', 'C':'2', 'T':'3', 'G':'4', 'a':'5', 'c':'6', 't':'7', 'g':'8'}
    b = {'1':'T', '2':'G', '3':'A', '4':'C', '5':'t', '6':'g', '7':'a', '8':'c'}
    reverse = [a[x] for x in seq][::-1]
    complement = ''.join([b[x] if x in b else x for x in reverse])
    return complement

  def to_fasta(self, genome):
    self.__try_load()
    data = []
    with open(genome) as fr:
      for l in fr.readlines():
        if l.startswith('>'):
          data.append([l.strip()[1:].split(' ')[0]])
        else:
          data[-1].append(l.strip())
    seqs = dict([(x[0], ''.join(x[1:])) for x in data])
  
    def store(file, feature):
      with open(file, 'w') as fw:
        for c, _, f, i, e, _, s, _, a in self.valid:
          if f == feature:
            id = dict([x.split('=') for x in a.split(';')])['ID']
            fw.write(f'>{id}\n{self.__strand(seqs[c][int(i)-1:int(e)], s)}\n')

    store(self.name+'_genes.fna', 'gene')
    store(self.name+'_mrnas.fna', 'mRNA')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: curl -s https://raw.githubusercontent.com/MiqueiasFernandes/bioinformatics/master/gff.py | python3 - genes.gff [newfile] [genome.fasta] [ptnas.txt]")
    else:
      file = sys.argv[1]
      gff = GFF(file)
      gff.ptnas = sys.argv[4] if len(sys.argv) > 4 else None
      print('Loading gff', file)
      gff.stats()
      if len(sys.argv) > 2:
        print('Persisting gff')
        gff.persist(sys.argv[2])
      if len(sys.argv) > 3:
       print('Storing fasta of genes and mRNAs')
       gff.to_fasta(sys.argv[3])
      print('all done.')
