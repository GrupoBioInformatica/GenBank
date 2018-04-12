from Bio import pairwise2 as pw2

first_seq = 'ATCGGATACCAA'
second_seq = 'ATCGGATTACCAA'

global_align = pw2.align.globalxx(first_seq, second_seq) #https://www.biostars.org/p/208540/
format_align = pw2.format_alignment(*global_align[0]) #help(pw2)

print(global_align[0])
print(format_align)

seq_length = min(len(first_seq), len(second_seq))
matches = global_align[0][2]
percent_match = (matches / seq_length) * 100

print("Percentage match: ", percent_match)
