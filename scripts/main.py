import random

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

def Mutate_DNA(seq): #פונקציה שמקבלת רצף דנא ומחליפה נוקלאוטיד לנוקלאוטיד אחר במקום רנדומלי
  new_seq = ''
  base_list = ['A', 'C', 'G', 'T']
  seq = seq.upper()
  ran_place = random.randrange(0, len(seq)) #הגרלת מקום רנדומלי
  ran_nuc = seq[ran_place] #הנוקלאוטיד במקום הרנדומלי שנבחר
  base_list.remove(ran_nuc) #מחיקת הנוקלאוטיד שנבחר מרשימת הנוקלאוטיד
  ran_base = random.randrange(0,len(base_list)) #הגרלת בסיס רנדומלי
  new_nuc = base_list[ran_base] 

  seq_1 = seq[0 : ran_place]
  seq_2 = seq[ran_place + 1 :] 
  new_seq = seq_1 + new_nuc + seq_2 #יצירת הרצף עם הנוקלאוטיד החדש
  return new_seq

def Insert_DNA(seq):
    base_list = ['A', 'C', 'G', 'T']
    seq = seq.upper()
    ran_place = random.randrange(0, len(seq))  # הגרלת מקום רנדומלי
    ran_nuc = seq[ran_place] #הנוקלאוטיד במקום הרנדומלי שנבחר
    ran_base = random.randrange(0, len(base_list))
    new_nuc = base_list[ran_base]

    seq_1 = seq[0: ran_place]
    seq_2 = seq[ran_place :]
    new_seq = seq_1 + new_nuc + seq_2
    return new_seq

def Delete_DNA(seq):
    seq = seq.upper()
    ran_place = random.randrange(0, len(seq))  # הגרלת מקום רנדומלי
    seq_1 = seq[0: ran_place]
    seq_2 = seq[ran_place + 1:]
    new_seq = seq_1 + seq_2
    return new_seq

RNA_codon_table = {}
codon_file = open('data/codon_AA.txt', 'r')
seq_file = open('data/p53_sequence.fa', 'r')
OG_dna = ''
Read_dict(codon_file)   

for line in seq_file:
  if line[0] != ">":
    line = line.rstrip('\n')
    OG_dna += line

remain = len(OG_dna) % 3 
OG_dna = OG_dna[0:len(OG_dna) - remain]
OG_rna = DNA_RNA_Cod(OG_dna)
OG_prot = RNA_prot(OG_rna)

mutated_dna = OG_dna

for i in range(0, 3): #יוצר מוטציה נקודתית 
  rand_mutation = random.randrange(0, 3)
  if rand_mutation == 0:
    mutated_dna = Mutate_DNA(mutated_dna)
  elif rand_mutation == 1:
    mutated_dna = Insert_DNA(mutated_dna)
  else:
    mutated_dna = Delete_DNA(mutated_dna)

remain = len(mutated_dna) % 3
mutated_dna = mutated_dna[0:len(mutated_dna) - remain]
mutated_rna = DNA_RNA_Cod(mutated_dna)
mutated_prot = RNA_prot(mutated_rna)
