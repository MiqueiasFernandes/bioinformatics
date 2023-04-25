import sys

try: from Bio import SeqIO
except: 
  print('Required: pip install biopython')
  exit()

def clean(file):
  seqs = SeqIO.parse(file, 'fasta')
  ok = []
  total = 0
  for seq in seqs:
    total += 1
    if seq.seq.endswith('*') or seq.seq.endswith('.') :
      ok.append(seq[:-1])
    elif str(seq.seq).isalpha():
      ok.append(seq)
    else:
      print('ERRO:', seq.id, '==> ', ','.join([x for x in seq.seq if not x.isalpha()])[:5])
  print('Total:', total)
  if all(['[protein_id=' in x.description for x in ok]):
    ver = []
    for x in ok:
      x.id = x.description.split("[protein_id=")[1].split("]")[0]
      if x.id in ver:
        raise Exception(f'ERROR: Duplicated `protein_id` {x.id}!')
      else:
        ver.append(x.id)
    with open(file+'_protein_id.txt', 'w') as f:
      f.writelines([f'{x}\n' for x in ver])
  final = SeqIO.write(ok, file+'_clean.faa', 'fasta')
  print('Clean:', final)
  print('Wrong:', total-final)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: curl -s https://raw.githubusercontent.com/MiqueiasFernandes/bioinformatics/master/ptnaClean.py | python3 - ptnas.faa")
    else:
      clean(sys.argv[1])
      