# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:54:58 2023

@author: Rabia KAŞIKCI
"""

import pandas as pd
import numpy as np

codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
        'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W',
}
aminoacids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y', 'V', 'X']

nuc_list=["A","T","G","C"]

def find_aminoacid(list_aa):
    aa_list=[]

    for i in range(len(list_aa)):
        for j in range(len(list_aa[i])):
            #print(list_aa[i][j])
            aa_list.append(codon_table[list_aa[i][j]])
            
    return aa_list
    


def generate_sequences(matrix):
    if not matrix:
        return [""]
    
    sequences = []
    current_element = matrix[0]

    for item in current_element if isinstance(current_element, list) else [current_element]:
        sequences.extend([item + seq for seq in generate_sequences(matrix[1:])])

    return sequences

new_codon = [[] for _ in range(3)]  

def generate_mutation(codon,mutation_num):
    new_list=[]
    if mutation_num==1:
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            new_codon[i] = [eleman for eleman in nuc_list if eleman != codon[i]]
            new_codon[j] = codon[j]
            new_codon[k] = codon[k]
            a=generate_sequences(new_codon)
            new_list.append(a)
            
    if mutation_num==2:
        for i in range(3):
            j = (i + 1) % 3  # Bu, döngüyü sağlamak için kullanılır
            k = (i + 2) % 3
            new_codon[i] = [eleman for eleman in nuc_list if eleman != codon[i]]
            new_codon[j] = [eleman for eleman in nuc_list if eleman != codon[j]]
            new_codon[k] = codon[k]
            a=generate_sequences(new_codon)
            new_list.append(a)
            
    if mutation_num==3:
        for i in range(1):
            new_codon[i] = [eleman for eleman in nuc_list if eleman != codon[i]]
            new_codon[i+1] = [eleman for eleman in nuc_list if eleman != codon[i+1]]
            new_codon[i+2] = [eleman for eleman in nuc_list if eleman != codon[i+2]]          
            a=generate_sequences(new_codon)
            new_list.append(a)
            
    return new_list




def generate_mut_table(codon,pro_list,mutation_matrix):
    x_index=aminoacids.index(codon_table[codon])
    for i in range(len(pro_list)):
        y_index=aminoacids.index(pro_list[i])
        mutation_matrix[x_index][y_index]=mutation_matrix[x_index][y_index]+1
        
    
    

aa_matrix1=np.zeros((21,21))
aa_matrix2=np.zeros((21,21))
aa_matrix3=np.zeros((21,21))


for i in nuc_list:
    for j in nuc_list:
        for k in nuc_list:
            codon=i+j+k
            #print(codon)     
            #1 mut

            pro_list_1=find_aminoacid(generate_mutation(codon,1))
            generate_mut_table(codon,pro_list_1,aa_matrix1)
            #2 mut
            pro_list_2=find_aminoacid(generate_mutation(codon,2))
            generate_mut_table(codon,pro_list_2,aa_matrix2)
            #3 mut           
            pro_list_3=find_aminoacid(generate_mutation(codon,3))
            generate_mut_table(codon,pro_list_3,aa_matrix3)


            
    

sum_of_mut=(aa_matrix1+aa_matrix2+aa_matrix3)

result=pd.DataFrame(aa_matrix1,columns=aminoacids, index=aminoacids)
result.to_excel("Single_Mut_Result.xlsx")            
            

result=pd.DataFrame(aa_matrix2,columns=aminoacids, index=aminoacids)
result.to_excel("Double_Mut_Result.xlsx")            
         

result=pd.DataFrame(aa_matrix3,columns=aminoacids, index=aminoacids)
result.to_excel("Trible_Mut_Result.xlsx")            
            

result=pd.DataFrame(sum_of_mut,columns=aminoacids, index=aminoacids)
result.to_excel("All_Mut_New_Result.xlsx")            
            
            

        