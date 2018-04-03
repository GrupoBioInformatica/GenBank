from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO #Librerias para acceder a las bases de datos


Entrez.email = "rau_95@hotmail.com" #Es necesario poner el email, aunque no es necesario hacerse cuenta.

with Entrez.esearch(db="nucleotide", term="Escherichia") as handle:
    record = Entrez.read(handle)

#Con esearch se crea, digamos, la consulta. Con read se ejecuta y devuelve datos de la bbdd.
#Parametros que he usado para esearch (hay mas):
    #db : la base de datos. Hay varias: nucleotide, gene, protein, pubmed, taxonomy,...
        #se pueden consultar, por ejemplo, en https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
    #term : terminos que busca. se puede indicar tambien el campo, ponendo valor[campo].
        #por ejemplo, "asthma AND DNA[TITL]" (buscaria asthma en cualquier campo y DNA en el titulo.
    

with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=record["IdList"][0]) as handle:
    seq_record = SeqIO.read(handle, "gb") #using "gb" as an alias for "genbank"

#Funciona de forma parecida a lo de arriba. efetch crea la consulra, SeqIO.read la ejecuta y recupera los datos.
#Parametros para efetch (hay mas):
    #db : la base de datos
    #rettype : el tipo de datos que se recuperan. (gb, fasta, xml, null,...)
    #retmode : como guardamos los datos recibidos. (text, xml, html,...) Consultar estos dos campos en https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

print("%s with %i features" % (seq_record.id, len(seq_record.features))) #id es la id del gen, features no se exactamente lo que es.
print(seq_record)
print("Sequence: ", seq_record.seq) #.seq es la secuencia de nucleotidos en si (por ejemplo, AATGCCTA)
#coding_dna = Seq(seq_record.seq, IUPAC.unambiguous_dna)

print("===================================================================")
print("===================================================================")


#Transcription
#Tenemos el ADN
template_dna = seq_record.seq.reverse_complement()


messenger_rna = seq_record.seq.transcribe()
print("Sequence: ", messenger_rna)



 
#Translation


print("===================================================================")
print("===================================================================")


translation_messenger_rna = messenger_rna.translate(to_stop=True)
print("Sequence: ", translation_messenger_rna)

print("===================================================================")
print("===================================================================")




#messenger_rna = "AUGCUAUUAUAA"

tamanio = len(messenger_rna)
codones_fin = ["UAA","UAG","UGA"]
i = 0
aux = ""
nlist = 0
list_genes = [ ]

if tamanio > 2:
	aux = messenger_rna[0]
	aux += messenger_rna[1]
	aux += messenger_rna[2]
	i = 3
#A continuacion vamos a separar el genoma en genes 
	while i < tamanio:
		
		gen = ""
		gen_util = False
		while aux != "AUG" and i < tamanio:
			gen += aux[0]

			var1 = aux[2]
			var0 = aux[1]
			aux = aux[1]
			aux += var1
			aux += messenger_rna[i]
			i+=1

		if i < tamanio:
			if len(gen) == 0:
				gen = aux
			print("Sequence: ", gen)
			print("Util: ", gen_util)

			gen_aux = [gen, gen_util]
			list_genes += gen_aux
			gen = aux
			gen_util = True
			nlist += 1
			

		print("===================================================================")
		print("===================================================================")
		
		while aux != codones_fin[0] and aux != codones_fin[1] and aux != codones_fin[2] and i < tamanio:
			aux = messenger_rna[i]
			i+=1
			aux += messenger_rna[i]
			i+=1
			aux += messenger_rna[i]
			i+=1
			gen += aux
		print("Sequence: ", gen)
		print("Util: ", gen_util)
		print("Aminoacidos: ", Seq(gen).translate(to_stop=True))

		print("===================================================================")
		print("===================================================================")
		gen_aux = [gen, gen_util]
		list_genes += gen_aux
		nlist += 1
				

	

#Este es otro ejemplo. No se porque, sale un gen en el que no se nombra el asma directamente.
#A lo mejor coge relaciones que no aparecen en los datos recuperados, habra que seguir investigando.

# Paginas web consultadas:
# http://biopython.org/DIST/docs/api/Bio.Entrez-module.html
# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
# https://www.ncbi.nlm.nih.gov/books/NBK3837/
# https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
# https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# https://www.ncbi.nlm.nih.gov
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc16
# https://www.ncbi.nlm.nih.gov/books/NBK49540/
