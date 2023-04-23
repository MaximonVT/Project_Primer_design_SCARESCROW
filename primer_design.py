#!/usr/bin/env python3 #important so the script can run!!!
"""
This script creates forward and reverse primers for targeted genes from the plant Arabidopsis thaliana. All you need to do is edit or adapt the cDNA.py file. Any gene inserted as a library type should work and the system will then give the ideal recommended primers / or if there are none, than the non-recommended primers with the respective melting temperature.
Author: Maximon
Created: 23.03.2023
Last Update: 15.04.23
"""
from cDNA import * #import a cDNA-Sequence

start_coder = "ATG"#define the start codon
end_coder = ["TAG", "TGA", "TAA"]#define the end codon

def reverse_complement(seq): #here we define the command reverse_complement to use it later in a simpler command line
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #complement the aminos
    return ''.join([complement[base] for base in reversed(seq)])#returns the complemented seq

def calculate_melting_temp(primer_seq):#here we define the calculating for the melting temperatur to keep it later easier
    num_A = primer_seq.count('A')#counting letters
    num_C = primer_seq.count('C')
    num_G = primer_seq.count('G')
    num_T = primer_seq.count('T')
    melting_temp =2*(num_A + num_T) + 4*(num_G + num_C)# equation of temperature calculation
    return melting_temp

def design_primers(genome):#here we define the code for the primers, from the genome -> later sequence, because we dont use the whole
    for genome_name, genome_seq in genome.items():#Genome is one dictionary data format with names and corresponding contents
        
        best_found = 0#If we found the best, it will be set as 1
        pos_starter = -1#can never be -1 so it should search unless he cant find one, than the system will give us a -1
        pos_ender = -1
        selected_for_primer = ""#“” means empty, which is used to initialize the variable
        selected_rev_primer = ""
        forward_temp = 0
        reverse_temp = 0
        for i in range(len(genome_seq)): # search all aequence # i starts from 0, 1, 2, 3, .....(loop)
            primer_seq = genome_seq[i:]
            if primer_seq.startswith(start_coder) and pos_starter == -1:# == means compare whether both sides are equal # only if we find any starter, and not any starter is found yet, we do the following and primer_seq[-3:] in end_coder:
                pos_starter = i # when finding the first start_coder, contimue looking for the ender

            if primer_seq[:3] in end_coder and pos_starter > -1: # :3 means, till the 3, i.e., 0,1,2 only if we find any ender and the starter is also already found
                pos_ender = i+3 # e.g., if the first letter of end is at pos 90, the end position shall be 92, pos_ender = 93
                
                # if the primer is found do the following
                one_sequence = genome_seq[pos_starter:pos_ender]
                selected_for_primer = ""
                selected_rev_primer = ""
                
                for l in range(19, 24): # the possiblity of length: 19, 20, ..., 23
                    forward_primer = one_sequence[:l]
                    reverse_primer = one_sequence[-l:]
                    #reverse_primer_mirror = reverse_primer[::-1] # mirror the reverse primer
                    # the function of temperature calculation is called for both primers
                    forward_temp = calculate_melting_temp(forward_primer)
                    reverse_temp = calculate_melting_temp(reverse_primer)#_mirror)

                    # check the temperature of both primers
                    if forward_temp >= 55 and forward_temp <= 62:
                        selected_for_primer = forward_primer

                    if reverse_temp >= 55 and reverse_temp <= 62:
                        selected_rev_primer = reverse_primer#_mirror

                    # check if we have found the primers we want
                    if selected_for_primer == "" or selected_rev_primer == "": # if any primer is not found, continue looking 
                        continue
                    elif not abs(forward_temp - reverse_temp) <= 4:# if the temperature difference is not 4, continue looking #Elif means else if #Not means, when the following condition is not met #elif not a>?, means if a<=? #abs means absolute value
                        best_found = 0
                        continue
                    else:
                        best_found = 1 # if the primers are found and the temperature difference is also 4, stop looking and enable "best_found" for later use
                        break
                
                if not best_found: # if the primers are not found reset_all and continue looking
                    pos_starter = -1
                    pos_ender = -1

            if best_found:
                break # if the primers are found and the temperature difference is also 4, stop looking

        
        #print the results
        if selected_for_primer == "" or selected_rev_primer == "":
            print("There are no possible primers for genome: ", genome_name)
        else:
            reverse_primer_complement = reverse_complement(selected_rev_primer)
            print("For genome:", genome_name)
            print("\t forward primer: ", selected_for_primer, "temperature: ", forward_temp)
            print("\t reverse primer (mirror and complement): ", reverse_primer_complement, "temperature: ", reverse_temp)
            if best_found:
                print("\t the combination of those primers are recommended!")
            else:
                print("\t the combination of those primers are not recommended!")

# the code starts here
design_primers(Genome)