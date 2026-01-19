def DNA_RNA_Cod(seq): #פונקציה שמקבלת רצף מקודד והופכת אותו לרנא
  RNA_seq = ''
  seq = seq.upper()
  for nuc in seq: # לולאה שמחליפה את נוקלאוטיד T לנוקלאוטיד U
    if nuc == "T":
      RNA_seq += 'U'
    else:
      RNA_seq += nuc  
  return RNA_seq

def Read_dict(file): #פונקציה שקוראת את הקובץ של הקודונים וחומצות האמינו והופפכת אותו למילון
  global RNA_codon_table
  for line in file:
    line = line.rstrip('\n')
    (Codon, sep, Amino_Acid) = line.partition('\t') #הרדה של הקודון מהחומצת אמינו בקובץ
    #print(Amino_Acid)
    RNA_codon_table[Codon] = Amino_Acid

def RNA_prot(seq): #פוקציה שמקבלת רצף רנא ומתרגמת אותו לרצף חלבונים
  prot_seq = ''
  for i in range(0, len(seq), 3): #לולאה שרצה על הרצף בקפיצות של 3
    codon = seq[i : i+3] #כל קודון שווה לשלישייה מהרצף
    Amino = RNA_codon_table[codon] #מציאת החומצת אמינו באמצעות המילון
    prot_seq += Amino
  return prot_seq



RNA_codon_table = {}
codon_file = open('data/codon_AA.txt', 'r')
seq_file = open('data/p53_sequence.fa', 'r')