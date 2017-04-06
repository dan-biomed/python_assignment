### provided function 1
def get_sequences_from_file(fasta_fn):
    sequence_data_dict = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        description = record.description.split()
        species_name = description[1] + " " + description[2]
        sequence_data_dict[species_name] = record.seq
    return(sequence_data_dict)



###function 2: adapted from YouTube tutorial "Translation from DNA to Protein. By Hong Qin"
def translate_function(string_nucleotides): 
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] 
    aa_seq_string =""
    for i in range(0, len(string_nucleotides), 3):  
        codon = string_nucleotides[i:i+3]   
        aminoacid = mito_table.forward_table[codon] 
        aa_seq_string += aminoacid        
    return(aa_seq_string)




###funtion 3: alternative function
from Bio.Seq import Seq
def alt_translate(alt_aa_seq):
###Define your sequence to be translated
    coding_dna = Seq(alt_aa_seq)
    alt_aa_seq = coding_dna.translate(table="Vertebrate Mitochondrial", to_stop=True)
    return alt_aa_seq



###function 4: aa charge and ratio (3 functions)

def charged(aa_seq): 
	charged = ['R','K','D','E'] 
	count = 0 
    for aa in aa_seq: 
        if aa in charged:
            count += 1 
    return(count/len(aa_seq)) 



def polar(aa_seq): 
	polar = ['Q','N','H','S','T','Y','C','M','W']
	count = 0 
    for aa in aa_seq: 
        if aa in polar:
            count += 1 
    return(count/len(aa_seq)) 



def hydrophobic(aa_seq): 
	hydrophobic = ['A','I','L','F','V','P','G']
	count = 0 
    for aa in aa_seq: 
        if aa in hydrophobic:
            count += 1 
    return(count/len(aa_seq)) 










## MAIN

###this is to create the dictionary of parse extracted names and sequences
cytb_seqs = get_sequences_from_file("bears_cytb.fasta") 
bear_df = pd.read_csv("bears_data.csv") 
species_list = list(bear_df.species)

bear_df



### resetting the counter
row_counter = 0
### for each key and value (species name and seq) in the cytochrome B dictionary, use my alternative translate function to get proportion of aa type (charged, polar, or hydrophobic). The function needs to run on a string
for key,value in cytb_seqs.items(): 
    aa_seq = alt_translate(str(value)) 
    charged_aa = charged(aa_seq)
    polar_aa = polar(aa_seq)
    hydrophobic_aa = hydrophobic(aa_seq)
    ### set the value for each proportion in the dataframe. The for loop iterates over rows one at a time
    bear_df.set_value(row_counter, 'charged', charged_aa) 
    bear_df.set_value(row_counter, 'polar', polar_aa) 
    bear_df.set_value(row_counter, 'hydrophobic', hydrophobic_aa)
    row_counter = row_counter + 1 
bear_df



### import the matplotlib magic function. Apparently there are so many magic functions!
%matplotlib inline
import seaborn as sea
import matplotlib.pyplot as plt
sea.barplot(x='species', y='mass', data=bear_df)
plt.xticks(rotation=45)
#to rotate labels so I can actually see: http://stackoverflow.com/questions/26540035/rotate-label-text-in-seaborn-factorplot

### The biggest bear was the U. spelaeus, which was over 500 kgs average. The composition of aa types are not very different across species. Although, as a biologist, I'd think that there could be real differences in the hydrophobic aa between the highest and lowest proportions, however small. I didn't expect that there is 1.6% difference between H. malayanus and U. arctos, but I suspect that if we compared across more genes this differences would be less.
### What else is interesting about this bear? This seems like an open-ended question, so I'll just say that apparently it was discovered in Italy that Neanderthals used to worship these cave bears.
