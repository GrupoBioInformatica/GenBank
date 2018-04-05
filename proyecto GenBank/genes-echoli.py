from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO #Librerias para acceder a las bases de datos


Entrez.email = "rau_95@hotmail.com" #Es necesario poner el email, aunque no es necesario hacerse cuenta.

with Entrez.esearch(db="nucleotide", term="escherichia[orgn] AND complete genome[title]") as handle:
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
a = input("pausa organismo")
print("Sequence (DNA): ", seq_record.seq) #.seq es la secuencia de nucleotidos en si (por ejemplo, AATGCCTA)
#coding_dna = Seq(seq_record.seq, IUPAC.unambiguous_dna)

print("===================================================================")
print("===================================================================")

#Transcription
#Tenemos el ADN
template_dna = seq_record.seq.reverse_complement()


messenger_rna = seq_record.seq.transcribe()
print("Sequence (transcrito): ", messenger_rna)


#Translation


print("===================================================================")
print("===================================================================")

translation_messenger_rna = messenger_rna.translate(to_stop=True)
print("Sequence (traducido): ", translation_messenger_rna)

print("===================================================================")
print("===================================================================")




#messenger_rna = "AUGCUAUUAUAA"

tamanio = len(messenger_rna)
a = input("pausa")

codones_fin = ["UAA","UAG","UGA"]
i = 0
aux = ""
nlist = 0
list_genes = [ ]

if tamanio > 2:
    aux = str(messenger_rna[i:i+3])
    #A continuacion vamos a separar el genoma en genes 
    while i < tamanio-3:

        gen = ""
        gen_util = False

        if aux != "AUG":

            while aux != "AUG" and i < tamanio-3:
                gen += aux[0]

                i += 1
                aux = str(messenger_rna[i:i+3]) #Queremos mirar de uno en uno?
            print("===================================================================")
            print("===================================================================")


            print("Sequence: ", gen)
            print("Util: ", gen_util)

            print("===================================================================")
            print("===================================================================")

        elif aux == "AUG":
            gen_util = True
            while aux != codones_fin[0] and aux != codones_fin[1] and aux != codones_fin[2] and i < tamanio-3:
                gen += aux

                i += 3
                aux = str(messenger_rna[i:i+3])
            gen += aux
            print("===================================================================")
            print("===================================================================")


            print("Sequence: ", gen)
            print("Util: ", gen_util)
            print("Aminoacidos: ", Seq(gen).translate(to_stop=True))

            print("===================================================================")
            print("===================================================================")


        list_genes += [gen, gen_util]
        nlist += 1

print("Numero de genes: ", nlist)
print("Tamanio: ", tamanio)
print("i: ", i)


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
