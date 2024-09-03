'''
Hiruni Anjana
Final Project Assignment
NCBI BLAST parser
2023/03/10
'''

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt

class BlastParser:

    def __init__(self, sequence):
        self.sequence = sequence

    @staticmethod
    #Method to output the type of the FASTA sequence
    def get_Seq_Type(sequence1):
        #List of amin acids
        aa = ["R", "N", "D", "B", "E", "Q", "Z", "H", "I", "L", "K", "M", "F", "P", "S", "W", "Y", "V"]
        #Read each base of the sequence
        for base in sequence1:
            #Check whether the sequence is mRNA, Amino acid or DNA
            if base not in aa and "U" in sequence1:
                return 'mRNA'
            elif base in aa:
                return 'Amino acid'
            elif base not in aa and "U" not in sequence1:
                return 'DNA'

    #Method to count the number of sequences for each type
    @staticmethod
    def count_sequences(fasta_file):
        #Initialize the counters
        dna_count = 0
        rna_count = 0
        aa_count = 0

        #Read the fasta file
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence = str(record.seq)

                #Check the type of sequence
                if set(sequence.upper()).issubset({'A', 'C', 'G', 'T'}):
                    dna_count += 1
                elif set(sequence.upper()).issubset({'A', 'C', 'G', 'U'}):
                    rna_count += 1
                elif set(sequence.upper()).issubset(
                    {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}):
                    aa_count += 1

        #Return the counts as a dictionary
        return {'DNA': dna_count, 'RNA': rna_count, 'Amino acid': aa_count}

    #A method to run BLAST searches on several sequence inputs in a single FASTA file
    def run_blast(self):
        #Read the FASTA file
        with open("fasta_file.fasta", "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                #Check the sequence type
                sequence_type = self.get_Seq_Type(record)

            #Run an appropriate BLAST algorithm type
            if sequence_type == 'DNA':
                result_handle = NCBIWWW.qblast('blastn', 'nt', record.seq)
            elif sequence_type == 'RNA':
                result_handle = NCBIWWW.qblast('blastn', 'nt', record.seq)
            else:
                result_handle = NCBIWWW.qblast('blastp', 'nr', record.seq)

            #Save the BLAST output for each sequence in XML format
            with open('record.id.xml', 'w') as out_handle:
                out_handle.write(result_handle.read())

    #Method to search for a given sequence motif in different BLAST hit sequences
    def search_motif(self,xml_file, motif):
        #Take the different BLAST records
        blast_records = NCBIXML.parse(open(xml_file))
        #Take attributes of each BLAST hit
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    #Check for the given sequence motif
                    if motif in hsp.sbjct:
                        print("Hit ID:", alignment.hit_id)
                        print("HSP sequence:", hsp.sbjct)

    #Method to parse through the BLAST hits from an XML file and output the sequence details, length, E-value, and the alignment
    def parse_blast(self, xml_file):

        from Bio.Blast import NCBIXML
        #Open the xml file
        with open(xml_file, 'r') as output:
            blast_record = NCBIXML.read(output)
            #Open a new text file to save the output
            output_file = open("Results.txt", "w")

            #Write the attributes of each blast hit to a text file
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    output_file.write("****Alignment****" )
                    output_file.write("\nAlignment length: \n" + str(alignment.length))
                    output_file.write("\nE-value: \n" + str(hsp.expect))
                    output_file.write("\nAlignment: \n" + str(alignment))
        #Close the file
        output_file.close()

    #Method to show the E-value distribution
    def evalue_histogram(self,xml_file):
        evalues = []
        #Open the xml file
        result_handle = open(xml_file)
        blast_records = NCBIXML.parse(result_handle)
        #Read each blast records
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    evalues.append(hsp.expect)
        result_handle.close()
        #Plot the histogram of E-values
        plt.hist(evalues, bins=50)
        plt.title("E-value Distribution")
        plt.xlabel("E-value")
        plt.ylabel("Frequency")
        plt.show()

    #Method that contains a user-given E-value threshold to filter out the BLAST hits in a BLAST output XML file
    def filter_blast_output(self,xml_file, evalue_threshold):
        from Bio.Blast import NCBIXML
        #Open the xml file
        with open(xml_file, 'r') as output:
            blast_record = NCBIXML.read(output)
            #Open a new text file to save the output
            output_file1 = open("E_value.txt", "w")

            #Write the attributes of each blast hit to a text file
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    #Check whether the BLAST hits less than the given E-value threshold
                    if hsp.expect < float(evalue_threshold):
                        output_file1.write("****Alignment****")
                        output_file1.write("\nBlast hit title: \n" + str(alignment.title))
                        output_file1.write("\nAlignment length: \n" + str(alignment.length))
                        output_file1.write("\nE-value: \n" + str(hsp.expect))
        #Close the file
        output_file1.close()















