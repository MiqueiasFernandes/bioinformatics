

import os, sys
import pymysql.cursors


	
def process(sql, arquivo, one, data):
	print('conectando ao banco de dados ...')
	try:
		connection = pymysql.connect(host='localhost',
									 user='root',
									 password='123',
									 db='guava',
									 charset='utf8mb4',
									 cursorclass=pymysql.cursors.DictCursor)
	except:
		print("ERRO: Impossivel conectar a base de dados.")
	else:
		result = None
		print("conectado com sucesso ... ")
		try:
			print("executando: " + sql + " " + str(data))
			with connection.cursor() as cursor:
				cursor.execute(sql, data)
				if one:
					result = cursor.fetchone()
				else:
					result = cursor.fetchall()
				#fo = open(arquivo, "w")
				if result:
					fo = open(arquivo, "w")
					print("salvando em: " + arquivo)
					if one:
						print("modo one")
						fo.write(str(result) + '\n')
					else:
						print("modo all")
						for res in result:
							fo.write(str(res) + '\n')
					#fo.close()
				else:
					print("o resultado não foi obtido!")
					fo = open(arquivo + '.err', "w")
					print("salvando ERRO em: " + arquivo + '.err')
					fo.write('o resultado não foi obtido!\n' + sql + '\n' + str(data))
				fo.close()
		except:
			print("ERRO: impossivel consultar tupla " + str(data) + ' query: ' + sql)
	print("terminado!")
	return(result)
	
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use: ./pysql sql arquivo all SQLarg1 SQLarg2 ... SQLargn")
    else:
        process(sys.argv[1], sys.argv[2], len(sys.argv) < 4 or not sys.argv[3] == 'all', sys.argv[4:])

