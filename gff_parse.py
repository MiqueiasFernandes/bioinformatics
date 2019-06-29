
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import hashlib

class Entry:
    def __init__(self, raw, igGen=(lambda e,d: e.name), data = None):
        self.raw = raw
        self.seqid = raw[0]
        self.sorce = raw[1]
        self.type = raw[2]
        self.start = int(raw[3])
        self.end = int(raw[4])
        self.score = raw[5]
        self.strand = raw[6]
        self.frame = raw[7]
        self.attributes = {k.split('=')[0]: k.split('=')[1] for k in raw[8].split(';')[:-1]}
        self.name = self.attributes['ID'] if 'ID' in self.attributes else (self.attributes['Name'] if 'Name' in self.attributes else '')
        self.parents = self.attributes['Parent'].split(',') if 'Parent' in self.attributes else []
        self.id = igGen(self, data)
        self.entryChilds = []
        self.entryParents = []
        self.size = 1 + self.end - self.start
        self.mark = False
        self.seq = None
        self.phase = 1 if self.frame == '1' else (2 if self.frame == '2' else 0)

    def clearMark(self):
        self.mark = False
        for c in self.entryChilds:
            c.mark = False
    
    def asGFF(self, alias=lambda e: "Alias=%s;" % e.attributes['Alias'] if 'Alias' in e.attributes else '', childs=True):
        self.mark = True
        return "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tID=%s;%s%s%s" % (
            self.seqid,
            self.sorce,
            self.type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.id,
            'Parent=%s;' % ','.join([e.id for e in self.entryParents]) if len(self.entryParents) > 0 else '',
            alias(self) if not alias is None else '',
            "\n" if childs and any([not e.mark for e in self.entryChilds]) else ''
        ) + ('\n'.join([e.asGFF(alias) for e in 
                       sorted(self.entryChilds, key=lambda e: 
                              1 if e.type == 'mRNA' else 
                              ( 2 if e.type == 'five_prime_UTR' else 
                               ( 3 if e.type == 'exon' else 
                                ( 4 if e.type == 'intron' else 
                                 ( 5 if e.type == 'CDS' else 6))))) if not e.mark]) if childs else '')
    
    def asFeature(self, seqid):
        return {'name': self.id, 
                'source': self.source, 
                'start': self.start, 
                'end': self.end, 
                'score': self.score,  
                'strand' : self.strand, 
                'frame' : self.frame,
                'sequence_id': seqid}
    
    def asSequence(self, fasta):
        if seq in None:
            raise 'Erro seq not imported, call writeFasta first'
        return {'name': self.id, 
                'size': self.size, 
                'is_contig': not 'N' in self.seq.upper(), 
                'type': 'GENE' if self.type == 'gene' else (
                    'MRNA' if self.type == 'mRNA' else (
                        'CDS' if self.type == 'CDS' else (
                            'EXON' if self.type == 'exon' else (
                                'INTRON' if self.type == 'intron' else (
                                    'FIVE_UTR' if self.type == 'five_prime_UTR' else (
                                        'THREE_UTR' if self.type == 'three_prime_UTR' else 'NUCL')))))),
                'sequence': self.seq
               }
    
    def hasType(self, t, childs=True):
        return self.type == t or any([ e.hasType(t, childs) for e in (self.entryChilds if childs else self.entryParents) ])
    
    def countType(self, t):
        return (1 if self.type == t else 0) + sum([e.countType(t) for e in self.entryChilds])
    
    def sizeOfType(self, t):
        return (self.size if self.type == t else 0) + sum([e.sizeOfType(t) for e in self.entryChilds])
    
    def writeGFF(self, file):
        file.write(self.asGFF(alias=None) + '\n')
    
    def setSequence(self, fasta=None):
        if self.seq is None and not fasta is None:
            s = fasta[self.seqid][self.start-1:self.end].seq
            if self.strand == '-':
                s = s.reverse_complement()
            self.seq = str(s)
        return self.seq
        
    def writeFasta(self, file=None, cols=80):
        if not file is None:
            file.write('>%s|%s|%d:%d\n%s\n' % ( 
                self.id, self.seqid, self.start, self.end, 
                '\n'.join([ self.seq[i:i+cols] for i in range(0, len(self.seq), cols) ])
                      ))
            
    def getMD5(self):
        return hashlib.md5(self.seq.encode('utf-8')).hexdigest()
        
    def __str__(self,):
        return self.asGFF(childs=False)
    def __repr__(self,):
        return str(self)


    

class GFF:
    def __init__(self, file, fasta, reload=True, create=True):
        if file is None:
            return
        print('loading %s ... ' % file)
        self.file_gff = file
        self.file_fasta = fasta
        ids = {}
        def idGenerator(entry, data):
            if not entry.type in data:
                data[entry.type] = 1
            else:
                data[entry.type] += 1
            return entry.type + str(data[entry.type])
        def orig(entry, data):
            return entry.name
        
        fn = idGenerator if create else orig

        self.gff = { k.id: 
        k for k in [ 
            Entry(l.strip().split('\t'), fn, ids) for l in set(open(file).readlines()) if l.count('\t') == 8 and not l.startswith('#')
                ] }
        print('%d entries importadas ...' % len(self.gff))

        self.genes = self.parseGenes()
        self.mrnas = self.parsemRNA()
        self.exons = self.parsemRNAResource('exon', 3)
        self.introns = self.parseIntrons()
        self.cdss = self.parsemRNAResource('CDS')
        self.futrs = self.parsemRNAResource('five_prime_UTR')
        self.tutrs = self.parsemRNAResource('three_prime_UTR')
        self.fasta = None

        genesClean = self.genes[0]
        geneRem = []
        mrnaRem = [x for x in self.mrnas[0] if len(x.entryChilds) < 1]
        mrnaClean = [x for x in self.mrnas[0] if len(x.entryChilds) > 0]
        for mrna in mrnaRem:
            g = mrna.entryParents[0]
            g.entryChilds.remove(mrna)
            if len(g.entryChilds) < 1:
                genesClean.remove(mrna.entryParents[0])
                geneRem.append(mrna.entryParents[0])
        if len(geneRem):
            print("genes removidos %d ..." % len(geneRem))
        if len(mrnaRem):
            print("mrnas removidos %d ..." % len(mrnaRem))
        if len(self.exons[2]):
            print("exons removidos %d ..." % len(self.exons[2]))            
        if len(self.introns[2]):
            print("introns removidos %d ..." % len(self.introns[2])) 
        if len(self.cdss[2]):
            print("cds removidos %d ..." % len(self.cdss[2]))
        if len(self.futrs[2]):
            print("5'UTR removidos %d ..." % len(self.futrs[2]))
        if len(self.tutrs[2]):
            print("3'UTR removidos %d ..." % len(self.tutrs[2]))

        self.gene = genesClean 
        self.mrna = mrnaClean 
        self.exon= self.exons[0] 
        self.intron= self.introns[0] 
        self.cds= self.cdss[0] 
        self.five_prime_UTR = self.futrs[0]
        self.three_prime_UTR = self.tutrs[0]
        self.all_entries = []
        self.all_entries.extend(self.gene)
        self.all_entries.extend(self.mrna)
        self.all_entries.extend(self.exon)
        self.all_entries.extend(self.intron)
        self.all_entries.extend(self.cds)
        self.all_entries.extend(self.five_prime_UTR)
        self.all_entries.extend(self.three_prime_UTR)
        self.dic = {e.id : e for e in self.all_entries}
        self.proteins = {}
        self.removed = {
                'gene':geneRem, 
                'mrna': mrnaRem, 
                'exon': self.exons[2], 
                'intron': self.introns[2], 
                'cds': self.cdss[2], 
                'five_prime_UTR': self.futrs[2],
                'three_prime_UTR': self.tutrs[2]
            }
            
        self.added= {
                'intron': self.introns[3]
            }

        if reload and (len(geneRem) + len(mrnaRem) + len(self.exons[2]) + len(self.introns[3]) + len(self.introns[2]) + len(self.cdss[2]) + len(self.futrs[2]) + len(self.tutrs[2])) > 0:
            new = file + '.new.gff3'
            self.genes = genesClean
            print("gerando novo gff alterado em %s ..." % new)
            self.storeGFF(new)
            print('*** RELOAD ***')
            self.__init__(new, fasta)
        elif not fasta is None:
            print('loading %s ...' % fasta)
            self.fasta = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
            print("set sequence for %d entries..." % len(self.all_entries))
            for e in self.all_entries:
                e.setSequence(self.fasta)
            self.proteins = self.getProteins()
        print('GFF %s loaded! ' % file)

    def parseGenes(self):
        entries = self.gff.values()
        cont = 1
        nGenes = []
        dic = {}
        for e in sorted(filter(lambda k: k.type == 'gene', entries), key=lambda x: (int(x.seqid.split('.')[2]) if '.' in x.seqid else x.seqid, x.start)):
            g = Entry([
                e.seqid, 
                e.sorce, 
                'gene', 
                str(e.start), 
                str(e.end), 
                e.score, 
                e.strand, 
                e.frame, 
                "ID=g%d;Alias=%s;" % (cont, e.name)])
            dic[e.name] = g
            nGenes.append(g)
            cont += 1
        print('%s genes parsed ...' % len(nGenes))
        return nGenes, dic

    def parsemRNA(self):
        entries = self.gff.values()
        genes = self. genes[1]
        cont = 1
        nmRNA = []
        dic = {}
        gnsc = {}
        for e in sorted(filter(lambda k: k.type == 'mRNA', entries), key=lambda x: (int(x.seqid.split('.')[2]) if '.' in x.seqid else x.seqid, x.start)):
            if len(e.parents) < 1:
                print('erro: %s (%s) nao tem parent!' % (e.id , e.name), file=sys.stderr)
            g = genes[e.parents[0]] if len(e.parents) > 0 and e.parents[0] in genes else None
            nid = 'mrna%da' % cont
            if not g is None and g.id in gnsc:
                gnsc[g.id][1] += 1
                nid = 'mrna%d%s' % (gnsc[g.id][0], 'abcdefghijklmnopqrstuvxwyz'[gnsc[g.id][1]])
            else:
                gnsc['-' if g is None else g.id] = [cont, 0]
                cont += 1
            m = Entry([
                e.seqid, 
                e.sorce, 
                'mRNA', 
                str(e.start), 
                str(e.end), 
                e.score, 
                e.strand, 
                e.frame, 
                "ID=%s;Parent=%s;Alias=%s;" % (nid, '' if g is None else g.id, e.name)])
            if not g is None:
                g.entryChilds.append(m)
                m.entryParents.append(g)
            dic[e.name] = m
            nmRNA.append(m)
        print('%s mRNA parsed ...' % len(nmRNA))
        return nmRNA, dic
        
    def parsemRNAResource(self, typeR, minSize=1):
        entries = self.gff.values()
        mrnas = self.mrnas[1]
        cont = 1
        rs = []
        dic = {}
        removidos = []
        for e in sorted(filter(lambda k: k.type == typeR, entries), key=lambda x: (int(x.seqid.split('.')[2])  if '.' in x.seqid else x.seqid, x.start)):
            if e.size < minSize:
                removidos.append(e)
                continue
            ms = [mrnas[p] for p in filter(lambda x: x in mrnas, e.parents)]
            e = Entry([
                e.seqid, 
                e.sorce, 
                typeR, 
                str(e.start), 
                str(e.end), 
                e.score, 
                e.strand, 
                e.frame, 
                "ID=%s%d;Parent=%s;Alias=%s;" % (typeR, cont, ','.join([m.id for m in ms]), e.name)])
            dic[e.name] = e
            for m in ms:
                m.entryChilds.append(e)
            e.entryParents = ms
            rs.append(e)
            cont += 1
        print('%d %s parsed ...' % (len(rs) + len(removidos), typeR))
        return rs, dic, removidos

    def parseIntrons(self):
        entries = self.gff.values()
        mrnas = self.mrnas[1]
        ret = self.parsemRNAResource('intron')
        introns = ret[0]
        dic = ret[1]
        skiped = ret[2]
        add = []
        for gene in set([m.entryParents[0] for m in mrnas.values() if m.hasType('gene', False) and m.entryParents[0].type == 'gene']):
            if gene.countType('intron') > 0:
                continue
            t = set()
            for l in [ [y for y in x.entryChilds if y.type == 'exon'] for x in gene.entryChilds]:
                t = t.union(l)
            t = sorted(t, key=lambda x: x.start)
            ps = {x: 0 for x in set(','.join([k.attributes['Parent'] for k in t]).split(','))}
            ms = {m.id: m for m in filter(lambda k: k.type == 'mRNA', gene.entryChilds)}
            ins = {}
            for e in t:
                for p in e.parents:
                    if ps[p] < 1:
                        ps[p] = e.end + 1
                    else:
                        k = "%d-%d" % (ps[p], e.start - 1)
                        if k in ins:
                            ins[k].append(p)
                        else:
                            ins[k] = [p]
                        ps[p] = e.end + 1
            for i in sorted(ins, key=lambda x: int(x.split('-')[0])):
                pos = i.split('-')
                n = Entry([
                    gene.raw[0], 
                    'script', 
                    'intron',
                    pos[0], 
                    pos[1], 
                    '.', 
                    gene.raw[6], 
                    '.',"ID=intron%d;Parent=%s;" % (len(introns) + 1, ','.join(ins[i]))])
                if n.size < 1:
                    skiped.append(n)
                    continue
                dic[n.id] = n
                introns.append(n)
                add.append(n)
                for p in ins[i]:
                    ms[p].entryChilds.append(n)
                    n.entryParents.append(ms[p])
        if len(add) > 0:
            print("introns adicionados %d ..." % len(add)) 
        return introns, dic, skiped, add
    
    def storeGFF(self, file=None, prefix=''):
        fn = prefix + ((self.file_gff + '.new.gff3') if file is None else file)
        with open(fn, 'w') as f:
            f.write('##gff-version 3\n')
            cont = 0
            for gene in self.gene:
                gene.clearMark()
                gene.writeGFF(f)
                cont += 1
        print('%d genes stored at %s' % (cont, fn))
    
    def storeFasta(self, prefix='', cols=80):
        filePrefix = prefix
        def export(dt, tp):
            print('salvando %s em %sgenome.%s.fna ...' % (tp, filePrefix, tp))
            with open('%sgenome.%s.fna' % (filePrefix, tp), 'w') as file:
                for x in dt:
                    x.writeFasta(file, cols)
        export(self.gene, 'gene')
        export(self.mrna, 'mrna')
        export(self.exon, 'exon')
        print('salvando cds em %sgenome.cds.fna ...' % filePrefix)
        with open('%sgenome.cds.fna'% filePrefix, 'w') as file:
            for mrna in sorted(filter(lambda m: m.hasType('CDS'), self.mrna), key=lambda x: (int(x.id[4:-1]), x.start, x.end)):
                cds = [c for c in sorted(mrna.entryChilds, key=lambda x: x.end if x.strand == '-' else x.start, reverse=mrna.strand == '-') if c.type == 'CDS']              
                seq = ''.join([e.seq for e in cds])
                file.write('>%s|%s| %s\n%s\n' % ( 
                    mrna.id, mrna.seqid, ' + '.join(['%d:%d' % (e.start, e.end) for e in cds]), 
                    '\n'.join([ seq[i:i+cols] for i in range(0, len(seq), cols) ])))
        print('salvando proteinas em %sgenome.protein.faa ...' % filePrefix)
        ptnas = self.proteins
        with open('%sgenome.protein.faa' % filePrefix, 'w') as file:
            for ptna in sorted(ptnas, key=lambda x: (int(x.split('|')[1][1:]), int(x.split('|')[0][4:-1]), x.split('|')[0][-1])):
                seq = ptnas[ptna]
                file.write('>%s\n%s\n' % (ptna, '\n'.join([ seq[i:i+cols] for i in range(0, len(seq), cols) ])))


    def getProteins(self, toStop=False, stopCodon='.', whithPhase=False):
        ret = {
            "%s|%s" % (m.id, m.entryParents[0].id): 
            str(Seq(
                ''.join([
                    c.seq[0 if c.strand == '-' else c.phase: (len(c.seq) - c.phase) if c.strand == '-' else len(c.seq)] if whithPhase else c.seq for c in sorted(m.entryChilds, key=lambda x: x.end if x.strand == '-' else x.start, reverse=m.strand == '-') if c.type == 'CDS'
                ]), generic_dna).translate(to_stop=toStop))
            .replace('*', stopCodon) 
            for m in filter(lambda m: m.sizeOfType('CDS') > 2, self.mrna)
        }
        print("proteins parsed %d ..." % len(ret))
        return ret
    
    def minimize(self, mrnaObg):
        mrna_obg = [l.strip() for l in set(open(mrnaObg).readlines()) if len(l) > 3]
        threshold = min([m.entryParents[0].size for m in gff.mrna if m.id in mrna_obg])
        genes = set([m.entryParents[0] for m in gff.mrna if m.size >= threshold ])
        print('reduced %.2f%% genes [at %dpb]... ' % ((100 - (len(genes)*100/len(self.gene))), threshold))
        s = GFF(None, None)
        s.gene = genes
        s.storeGFF(self.file_gff + '.tmp.gff3')
        g = GFF(self.file_gff + '.tmp.gff3', None)
        g.storeGFF(self.file_gff + '.min.gff3')
        return GFF(self.file_gff + '.min.gff3', self.file_fasta)
        
    def stats(self):
        print('calc stats ...')
        def calcTam(sizes, classes=[[1, 100], [100, 500], [500, 1000], [1000, 2000], [2000, 3000], [3000, 10**5]]):
            return {
                '1.min': min(sizes),
                '2.max': max(sizes),
                '3.mean': (sum(sizes) / len(sizes)),
                '4.qtd': len(sizes),
                '5.size': sum(sizes),
                '6.class': {
                    str(n[0]) + '.>=' + str(n[0]) + 'pb': len([x for x in sizes if x >= n[0] and x < n[1]]) for n in classes
                }
            }

        def numDeTypePerType(entries, t, classes=[[0, 1], [1, 2], [2, 5], [5, 10], [10, 15], [15, 10**5]]):
            qtds = [x.countType(t) for x in entries]
            ret = { 
                str(n[0]+1) + '.>=' + str(n[0]): len([x for x in qtds if x >= n[0] and x < n[1]]) for n in classes 
            }
            ret['0.mean'] = (sum(qtds) / len(qtds))
            return ret
        
        def toTab(dic, pad=''):
            s = ''
            for k in sorted(dic, key=lambda x: int(x.split('.')[0])):
                s += pad + k.split('.')[1] + '\t'
                if type(dic[k]) == dict:
                    s += '\n' + toTab(dic[k], pad + '    ')
                elif type(dic[k]) == float:
                     s += ("%.2f" % dic[k])
                else:
                    s += str(dic[k])
                s += '\n'
            return s
        dic = {
            '0.gene': {
                '0.size': calcTam([g.size for g in self.gene]),
                '1.num': {
                    '0.mrna': numDeTypePerType(self.gene, 'mRNA'),
                    '1.exon': numDeTypePerType(self.gene, 'exon'),
                    '2.intron': numDeTypePerType(self.gene, 'intron'),
                    '3.cds': numDeTypePerType(self.gene, 'CDS'),
                },
                '2.completo': calcTam([g.size for g in self.gene if g.hasType('five_prime_UTR') and g.hasType('three_prime_UTR')]),
                '3.incompleto': calcTam([g.size for g in self.gene if not g.hasType('five_prime_UTR') or not g.hasType('three_prime_UTR')]),
                '4.vazio': calcTam([g.size for g in self.gene if not g.hasType('five_prime_UTR') and not g.hasType('three_prime_UTR')]),
            },
            '1.mrna': {
                '0.size': calcTam([m.size for m in self.mrna]),
                '1.cds': calcTam([m.sizeOfType('CDS') for m in self.mrna]),
                '2.num': {
                    '0.exon': numDeTypePerType(self.mrna, 'exon'),
                    '1.intron': numDeTypePerType(self.mrna, 'intron'),
                    '2.cds': numDeTypePerType(self.mrna, 'CDS'),
                },
                '3.exon:intron': '1:2',
                '4.exon:cds': '1:3'
            },
            '2.exon': {
                '0.size': calcTam([e.size for e in self.exon]),
                '1.redundant': 0 #####add
            },
            '3.intron': {
                '0.size': calcTam([i.size for i in self.intron]),
                '1.redundant': 0 #####add
            },
            '4.cds': {
                '0.size': calcTam([c.size for c in self.cds]),
                '1.redundant': 0 #####add
            },
            '5.5\'UTR': {
                '0.size': calcTam([f.size for f in self.five_prime_UTR])
            },
            '6.3\'UTR': {
                '0.size': calcTam([t.size for t in self.three_prime_UTR])
            },
            '7.proteins': {
                '0.size': calcTam([len(self.proteins[p]) for p in self.proteins])
            }
        }
        return dic, toTab(dic)
    
    def storeStats(self, file=None, prefix=''):
        fn = prefix + ((self.file_gff + '.stats') if file is None else file)
        with open(fn, 'w') as f:
            f.write(self.stats()[1])
        print('STATS salvo em ' + fn)
        
    def writeMD5(self, file=None, prefix=''):
        fn = prefix + ((self.file_gff + '.MD5') if file is None else file)
        def wmd5(file, entry):
            file.write('%s\t%s\n' % (entry.id, entry.getMD5()))
            for child in entry.entryChilds:
                wmd5(file, child)

        with open(fn, 'w') as f:
            for gene in self.gene:
                wmd5(f, gene)
            for p in self.proteins:
                f.write('%s\t%s\n' % (p, hashlib.md5(self.proteins[p].encode('utf-8')).hexdigest()))
        print('MD5 salvo em ' + fn)
        
    def persistFile(self, prefix):
        self.storeFasta( prefix=prefix)
        self.storeGFF( prefix=prefix)
        self.storeStats(prefix=prefix)
        self.writeMD5(prefix=prefix)
    
    def persistMySQL(self):
        return True
