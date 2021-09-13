# mRNA-Morton-Project
Project to study mRNA secondary structures and conserved sites within plant photosynthesis genes. 
This is the procedure followed to analyze the MFE and variation of the third codon site for 8 species and 9 different genes. 

We started by verifying we have extracted the correct sequences associated with the genes and species of interest. For the purpose of our initial analysis we looked at 8 different species (Zea mays, Agave Americanum, Arabidopsis thaliana. Pisum sativum, Nicotiana tabacum, Aster altaicus, Magnolia macrophylla, Nymphaea alba) and 9 genes (atpA, atpB, psaA, psaB, psbA, psbB, psbC, psbD, rbcL). We had codon variation data for 7 genes and did not have codon variation data for atpA and psbB. 

The ViennaRNA package was used in this analysis (https://www.tbi.univie.ac.at/RNA/). We used this link to download and install ViennaRNA on the computer terminal: https://www.tbi.univie.ac.at/RNA/#download. To go through the tutorial and understand how to use ViennaRNA on your terminal, follow this: https://www.tbi.univie.ac.at/RNA/tutorial/. The package was installed in my home directory and the executable programs after following the tutorial were found in /Users/shailafye/Tutorial/Progs. 

After verifying we had the correct sequence and installing ViennaRNA, the program MFE_of_windows.py was run using the each gene sequence for the corresponding species. To execute this program, run script two times through and then switch to the terminal and in /Users/shailafye/Tutorial/Progs directory run both executables listed in MFE_of_windows.py. Make sure that the files saved match the correct directory. The files: sequence_windows.txt and all_random_seq.txt should be saved in /Users/shailafye/Tutorial/Progs directory with your ViennaRNA executable programs. Run the two commands on the command line: RNAfold < sequence_windows.txt > output_original.txt and RNAfold < all_random_seq.txt > random_output.txt. This will calculate MFE for each window. Then it gets saved as a file for each run: output_original.txt and random_output.txt. After running this in the terminal, re-run MFE_of_windows.py twice and then the new MFE data from ViennaRNA will be parsed in our written program. The MFE_of_windows.py program parses the MFE data generated by RNAfold executable on the terminal from the ViennaRNA package. Our program extracts the MFE and averages the MFE of the 100 randomized sequence MFE for each window which corresponds to (random_output.txt). This data is then saved as original_mfe.txt and random_avg.txt. These two files will be overwritten next time you run so make sure to save the data separately. 

For the purposes of our analysis, we saved the data in an Excel file for each gene for each species. From here, we graphed the MFE for each window and compared it with the randomized and average MFE. We did not initially see much value or significance by comparing with the randomized and average MFE. 

After graphing all the genes MFE for the 8 species, along with the average of all 8 species, consensus sequences were also generated and run through the same procedure above for calculating the MFE. We wrote the program: consensus_final.py which generated a consensus sequence based on a specified parameter, which in this case was 0.5. Also, using ViennaRNA’s program RNALalifold, with the input of the alignment, a consensus sequence was also generated. The alignment for the 8 species for each gene was done using: https://www.ebi.ac.uk/Tools/msa/clustalo/. This alignment was also the input for our program consensus_final.py. After generating consensus sequences, the sequence was run through the window generator and MFE program on ViennaRNA.

After finishing the MFE calculations, we moved on to looking at the third codon site and it’s variation. To do this, a file containing the codons and corresponding variation number for 169 previously studied species was used. This file was previously generated by Professor Morton and also contained whether the site was a 2 or 4 fold degenerate site. This data was processed by Codon_Variation_with_GC_content.py to find the variation for each window size we had previously looked at in MFE_of_windows.py. For the program Codon_Variation_with_GC_content.py, the input is a previously generated excel file with codons and the variation for each 2 and 4 fold degenerate listed. Make sure this Excel file has no empty rows at the bottom and any empty codon cells should be replaced with XXX. The program parses through the excel file and outputs the average variation for 2 and 4 fold variation. It also outputs the variation for sites that are surrounded by 0, 1, or 2 G’s or C’s to see if there is any correlation between stability and MFE. We compared the variation for 2 and 4 fold variation for each window to see if there was any correlation between variation and MFE. 

One last program of importance is extractcodonsequence.py, which takes in the same excel file and outputs a list of codons and the sequence in FASTA format. This was compared as a check for the above analysis using https://www.ebi.ac.uk/Tools/msa/clustalo/ with the sequences of the 8 species for that gene. Any modifications were made to the excel file to shift any codons or delete any. 






