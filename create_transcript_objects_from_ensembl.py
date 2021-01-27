###################
# REVISION HISTORY:
###################


##########
# IMPORTS:
##########

import glob    # global?            - for finding many files in a directory that match a given pattern
import re      # regular expression - for changing contents of strings easily


############
# FUNCTIONS: 
############

def process_output_line():
    pass

############
# MAIN CODE:
############

try:
    outfile = open("transcript_vars_objects_from_ensembl.py", "w")
except:
    print("output file error on open...halting process")
    exit()

#######

file_list = glob.glob("Homo_sapiens_*_sequence.fa") # creates a list of filenames matching the pattern in CURRENT directory


for file_name in file_list:

    current_gene_and_transcript = ""
    current_sequence            = ""
    current_exon_number         = 1
    addition_string             = ""
    exon_length_list            = []
    previous_header_line        = None

    try:
        infile = open(file_name, "r")
    except:
        print("input file error on open...halting process")
        exit()

    for line in infile:
        
        if  line == "\n": # if the entire line is a carriage return, skip line
            continue      # continue command brings us directly to the next loop iteration

        line = re.sub(r'[^a-zA-Z0-9-> ]',r'',line) # remove any char that isn't letter, number, dash or >
        print(line)


        if  line[0] == ">":

            """
               we're processing a header line, so something has "changed" (or we're on line 1 of the input file):
            """

            split_line_list = line.split() # split line spaces as delimiter to define list elements

            if  previous_header_line == None: # we have no previous_header_line yet so we must be on very first line of the file (transcript_1 of exon_1 ...we can ASSUME this is an exon header, and NEVER a UTR)
                current_gene_and_transcript = split_line_list[0]
                current_exon_number         = 1
                previous_header_line        = line
                continue                # move on to the next line in the input file (skip all remaining code in this loop iteration)

            else:
                build_new_line = current_gene_and_transcript # current_gene_and_transcript is in format geneName-transcriptName eg PHPT1-202

                if  "exon" in previous_header_line:
                    build_new_line += "_EXON" + str(current_exon_number)
                    addition_string += build_new_line + " + "
                    exon_length_list.append(len(current_sequence))
                if  "utr3" in previous_header_line:
                    build_new_line += "_UTR_3"
                    addition_string += build_new_line + " + "
                    exon_length_list.append(len(current_sequence))
                if  "utr5" in previous_header_line:
                    build_new_line += "_UTR_5"
                    addition_string = build_new_line + " + " + addition_string
                    exon_length_list.insert(0, len(current_sequence))

                build_new_line += " = " + "'" + current_sequence + "'"

                build_new_line = re.sub(r'[-]',r'_',build_new_line) # change all - to _
                build_new_line = re.sub(r'[>]',r'',build_new_line)  # remove all >'s

                outfile.write(build_new_line)
                outfile.write("\n")
                

                current_sequence            = ""
                
                if  split_line_list[0] != current_gene_and_transcript: # we've moved on to a NEW transcript for the gene
                    outfile.write("\n")
                    addition_string = current_gene_and_transcript + "_FULL_SEQUENCE = " + addition_string
                    addition_string = re.sub(r'[-]',r'_',addition_string) # change all - to _
                    addition_string = re.sub(r'[>]',r'',addition_string)  # remove all >'s
                    addition_string = addition_string.rstrip("+ ")

                    current_junctions = ""
                    for i in exon_length_list:
                        current_junctions += str(i) + ","
                    current_junctions = current_junctions.rstrip(",")
                    
                    class_df_string = current_gene_and_transcript + " = Transcript(" + current_gene_and_transcript + "_NAME," + current_gene_and_transcript + "_FULL_SEQUENCE, [" + current_junctions +"])"
                    class_df_string = re.sub(r'[-]',r'_',class_df_string) # change all - to _
                    class_df_string = re.sub(r'[>]',r'',class_df_string)  # remove all >'s

                    name_string = current_gene_and_transcript + "_NAME = " + "'" + current_gene_and_transcript + "'"
                    name_string = re.sub(r'[-]',r'_',name_string) # change all - to _
                    name_string = re.sub(r'[>]',r'',name_string)  # remove all >'s
                    outfile.write(name_string)
                    outfile.write("\n")
                    
                    outfile.write(addition_string)
                    outfile.write("\n")
                    outfile.write(class_df_string)
                    outfile.write("\n")
                    outfile.write("\n")
                    current_exon_number  = 1
                    addition_string      = ""
                    exon_length_list     = []
                else:                                                  # we've moved on to the next exon/UTR in the same transcript for the gene
                    current_exon_number += 1

                current_gene_and_transcript = split_line_list[0]
                previous_header_line        = line


        if  line[0] in ["A", "T", "C", "G"]: # first char isn't a '>' so we are simply adding more letters to the sequence
            current_sequence += line

####################
# last line in file: ( ********* this should go in a subroutine with the code above!! ********* )

##ESPECIALLY COS ADDITION LINE FOR FINAL LINE IS NOT SHOWING UP!!

    build_new_line = current_gene_and_transcript # current_gene_and_transcript is in format geneName-transcriptName eg PHPT1-202

    if  "exon" in previous_header_line:
        build_new_line += "_EXON" + str(current_exon_number)
        addition_string += build_new_line + " + "
        exon_length_list.append(len(current_sequence))
    if  "utr3" in previous_header_line:
        build_new_line += "_UTR_3"
        addition_string += build_new_line + " + "
        exon_length_list.append(len(current_sequence))
    if  "utr5" in previous_header_line:
        build_new_line += "_UTR_5"
        addition_string = build_new_line + " + " + addition_string
        exon_length_list.insert(0, len(current_sequence))

    build_new_line += " = " + "'" + current_sequence + "'" # building the last line from current file here!

    build_new_line = re.sub(r'[-]',r'_',build_new_line) # change all - to _
    build_new_line = re.sub(r'[>]',r'',build_new_line)  # remove all >'s

    outfile.write(build_new_line)
    outfile.write("\n")

    current_sequence            = ""
    
    if  True: # force TRUE because this is the last line in the file #split_line_list[0] != current_gene_and_transcript: # we've moved on to a NEW transcript for the gene
        outfile.write("\n")
        addition_string = current_gene_and_transcript + "_FULL_SEQUENCE = " + addition_string
        addition_string = re.sub(r'[-]',r'_',addition_string) # change all - to _
        addition_string = re.sub(r'[>]',r'',addition_string)  # remove all >'s
        addition_string = addition_string.rstrip("+ ")

        current_junctions = ""
        for i in exon_length_list:
            current_junctions += str(i) + ","
        current_junctions = current_junctions.rstrip(",")
        
        class_df_string = current_gene_and_transcript + " = Transcript(" + current_gene_and_transcript + "_NAME," + current_gene_and_transcript + "_FULL_SEQUENCE, [" + current_junctions +"])"
        class_df_string = re.sub(r'[-]',r'_',class_df_string) # change all - to _
        class_df_string = re.sub(r'[>]',r'',class_df_string)  # remove all >'s

        name_string = current_gene_and_transcript + "_NAME = " + "'" + current_gene_and_transcript + "'"
        name_string = re.sub(r'[-]',r'_',name_string) # change all - to _
        name_string = re.sub(r'[>]',r'',name_string)  # remove all >'s
        outfile.write(name_string)
        outfile.write("\n")
        
        outfile.write(addition_string)
        outfile.write("\n")
        outfile.write(class_df_string)
        outfile.write("\n")
        outfile.write("\n")
        current_exon_number  = 1
        addition_string      = ""
        exon_length_list     = []
    else:                                                  # we've moved on to the next exon/UTR in the same transcript for the gene
        current_exon_number += 1
        outfile.write("addition_string")

    current_gene_and_transcript = split_line_list[0]
    previous_header_line        = line

# close current input file and loop back to next input file
    infile.close()

# processing is complete so close the only output file we opened
outfile.close()
