#################
# IMPORT SECTION:

import sqlite3
import re      # regular expression - for changing contents of strings easily


########################################
########## FILE and DB SECTION #########

"""
A sample experiments.txt is a tab-separated file follows
(do NOT include header - only include data and # symbols)

Gene, transcript_id and exon_no must exist on the GTF/DB


gene_name - transcript_id - exon_no

PHPT1	ENST00000371661.5	001
PHPT1	ENST00000371661.5	002
PHPT1	ENST00000371661.5	003
#
PHPT1	ENST00000497413.1	001
PHPT1	ENST00000497413.1	002
#

...NB - ALL transcripts MUST be followed by a line starting with #

The output file is a generated PYTHON file made up of variables and class instance declarations
using those variable declarations.
This generated Python file will be "included" (via exec command) into another Python program
"fastq_builder.py", where the Groups and Replicates can be chosen/defined based on
what is in the output_filename "fastq_transcript_vars_objects_from_SQL.py"
(and on what is in the non-SQL Transcript/Vars/Objects file generated from ensembl downloads)
"""


input_filename = "experiments.txt"
             
try:
    infile = open(input_filename, "r")
except:
    print("input file error on open...halting process")
    exit()

#
output_filename = "transcript_vars_objects_from_sql.py"

try:
    outfile = open(output_filename, "w")
except:
    print("output file error on open...halting process")
    exit()

#

conn = sqlite3.connect('fasta_gtf_v35.db')
c = conn.cursor()

########################################

experiments = []

for line in infile:

    if  line[0] != "#":                       # BUILD UP A LIST OF ENTRIES FOR THE EXPERIMENT
        split_string = line[:-1].split("\t")

        experiments.append(split_string)
    else:
        # new list/experiment                 # ONCE WE HAVE A COMPLETE LIST, PROCESS THE EXPERIMENT
##        print(experiments)

        gene_name = experiments[0][0]
        xcript_id = experiments[0][1]


        #######
        # NAME:
        #######
        exp_name = experiments[0][0] + "_" + experiments[0][1]
        for entry in experiments:
            exp_name += ("_" + entry[2])

##        print(exp_name)


        ################
        # FULL SEQUENCE:  (this will be a concatenation of all of the requested exons for this experiment):
        ################

        transcript_seq = ""   # will be combo of all exons for this experiment together
        list_of_exon_lens = []
        
        for entry in experiments:

            exon_no   = entry[2]

            # get To-From Coords and ChromosomeNo (Chrom)
            c.execute('SELECT CoordFrom, CoordTo, Chrom FROM EXONS_GTF WHERE GeneName = ? AND XcriptID = ? AND ExonNo = ?',(gene_name, xcript_id, exon_no))
            exons_gtf_rows = c.fetchone()
##            print("gene_name, xcript_id, exon_no", gene_name, xcript_id, exon_no)
##            print("\n exons_gtf_rows")
##            print(exons_gtf_rows)

            if  exons_gtf_rows == None: # row not found in EXONS_GTF, most likely meaning a <genome>.gtf definition has been updated
                                        # from previous definition of <genome>.gtf (transcript id, gene name or exon no)
                print("Fatal Error! GeneName + XcriptId + ExonNo Not Found in file EXONS_GTF! experiments.txt error: ",gene_name, xcript_id, exon_no)
                exit()                  # halt the execution of the program

            chrom_no = exons_gtf_rows[2]
            upper_index_from_K = str(int(exons_gtf_rows[0][:-3]) + 1)
            upper_index_to_K   = str(int(exons_gtf_rows[1][:-3]) + 1)
            #int_exon_no = int(exon_no)

##            print("fromK: ", upper_index_from_K)
##            print("toK:   ", upper_index_to_K)

            c.execute('SELECT Seq FROM FASTA_1000 WHERE Chrom = ? AND UpperIndex >= ? AND UpperIndex <= ?',(chrom_no, upper_index_from_K, upper_index_to_K))
            seq_parts = c.fetchall()
##            print("seq parts: ")
##            print(seq_parts)

            full_seq = ""
            for seq_part in seq_parts:
##                print("SEQ PART:  *************")
##                print(seq_part)
                full_seq += seq_part[0] # contents of seq_parts is a list of tuples, each tuple containing one string (we need to get to the string/s). This is a SQLITE "thing"
##            print("full_seq:")
##            print(full_seq)

            exon_span_K = int(upper_index_to_K) - int(upper_index_from_K) #  how many FASTA_1000 records does the exon data span over

            exon_range_from = exons_gtf_rows[0][-3:] # take only the last 3 chars
            exon_range_to   = exons_gtf_rows[1][-3:] # take only the last 3 chars
            exon_range_to   = str(exon_span_K) + exon_range_to

##            print("exon range: ")
##            print("from: ", exon_range_from)
##            print("to  : ", exon_range_to)

            exon_seq = full_seq[ (int(exon_range_from) - 1)  : int(exon_range_to)] # we need to start in 1 pos earlier because GTF range is base-1 but python is base-0
            
##            print("exon sequence: ", exon_seq)

            transcript_seq += exon_seq
            list_of_exon_lens.append(len(exon_seq))

        ####################
        # OBJECT DEFINITION:
        ####################

        # token definition (this is needed as I can't use . in Python variable names (they can still be in the value of course)
        token_exp_name = re.sub(r'[.]',r'_',exp_name) # change all - to _

        # name line
        name_string = token_exp_name + "_NAME = " + "'" + exp_name + "'"
##        print(name_string)

        # full seq line
        full_seq_line = token_exp_name + "_FULL_SEQUENCE = '" + transcript_seq + "'"
##        print(full_seq_line)

        # list of exon lengths
        current_junctions = ""
        for i in list_of_exon_lens:
            current_junctions += str(i) + ","
        current_junctions = current_junctions.rstrip(",")

        class_df_string = token_exp_name + " = Transcript(" + token_exp_name + "_NAME, " + token_exp_name + "_FULL_SEQUENCE, [" + current_junctions +"])"
        
##        print(class_df_string)

##        print("\n")

        outfile.write(name_string + "\n")
        outfile.write(full_seq_line + "\n")
        outfile.write(class_df_string + "\n")

        outfile.write("\n")
        
        experiments = [] # reset data for next loop

# output a space to make generated Python file more readable:
    

# close all flat files

infile.close()
outfile.close()

    

    
