from Project import BlastParser

if __name__ == "__main__":

    #Files
    fasta_file = "fasta_file.fasta"
    xml_file = "record.id.xml"
    motif = "ATGC"
    # file2 = "Results.txt"

    sequence1 = ""
    #Open the FASTA file
    with open("seq_q31.fasta") as file1:
        for line in file1:
            #Remove escape characters and FASTA header
            if line != "\n" and ">" not in line:
                #Remove the spaces and split the words with tabs
                line = line.strip()
                sequence1 += line

    #Calling the 1st method (Type of the FASTA sequence)
    type = BlastParser.get_Seq_Type(sequence1)
    print("Sequence type: ", type)

    #Calling the 2nd method (No. of sequences for each type)
    count = BlastParser.count_sequences(fasta_file)
    print("Number of sequences: ", count)

    #Creating object
    obj1 = BlastParser(fasta_file)
    #Calling 3rd method (Run the BLAST program)
    obj1.run_blast()

    #Calling the 4th method (Search for a given sequence motif)
    obj1.search_motif(xml_file, motif)

    #Calling the 5th method (Parse through BLAST hits)
    obj1.parse_blast(xml_file)

    #Calling the 6th method (E-value distribution)
    obj1.evalue_histogram(xml_file)

    #User input
    evalue_threshold = input("Enter the E-value threshold value: ")
    #Calling the 7th method (Filter out the BLAST hits)
    obj1.filter_blast_output(xml_file, evalue_threshold)

