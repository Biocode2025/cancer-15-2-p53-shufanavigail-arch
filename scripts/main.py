def DNA_RNA_Cod(seq): #פונקציה שמקבלת רצף מקודד והופכת אותו לרנא
  RNA_seq = ''
  seq = seq.upper()
  for nuc in seq: # לולאה שמחליפה את נוקלאוטיד T לנוקלאוטיד U
    if nuc == "T":
      RNA_seq += 'U'
    else:
      RNA_seq += nuc  
  return RNA_seq
