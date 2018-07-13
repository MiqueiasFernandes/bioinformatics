
import sys
import re
from subprocess import call

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("use: ./pysql input output")
	else:
		file = open(sys.argv[1], "r")
		incds = None
		inscf = None
		cont = 0
		errs = 0
		for line in file:
			if line.startswith(">"):
				incds = line.split(' ')[0][1:]
			elif line.startswith("pg.scf."):
				inscf = re.sub(r'\W*$', "", line)
			elif line.startswith("Sbjct"):
				if not (incds == None or inscf == None):
					spl = line.split(' ')
					try:
						s1 = spl[2]
						s2 = re.sub(r'\W*$', "", spl[4])
						ini = min(int(s1), int(s2))
						fim = max(int(s1), int(s2))
						#print("comp: " + incds + " / " + inscf + " / " + str(ini) + ":" + str(fim))
						
						call(["./scrit.sh", str(ini), str(fim), inscf, sys.argv[2] + str(cont + 1) + str(incds + '_' + inscf + '_' + str(ini) + '-' + str(fim))])
						cont += 1
					except:
						print("ERROR: ./scrit.sh" + " " + str(ini) + " " + str(fim) + " " + inscf + " " + str(incds + '_' + inscf + '_' + str(spl[2:]) + ' => '+ s1 + ':' + s2))
						errs += 1
				else:
					print("ERR: " + line[:-1] + incds + inscf)
					errs += 1
		#call(["./scrit.sh", '108666', '109206', 'pg.scf.563', 'saida.txt'])
		print("total " + str(cont) + " processadas ...")
		print("total " + str(errs) + " falharam ...")



