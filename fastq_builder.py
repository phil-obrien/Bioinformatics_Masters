####################
# CLASS DEFINITIONS:
####################

from random import randrange

class Transcript:

    def __init__(self, name, full_sequence, exon_len_list):
        self.name            = name
        self.full_sequence   = full_sequence
        self.exon_len_list   = exon_len_list

        self.junction_list = [] # this will be populated later in the code, derived from exon_len_list


####################################
# INCLUDE CODE FROM OUTSIDE SOURCES:
####################################

# The exec command brings in code from an external file. Used because "import" does not work for this functionality

try:
    exec(open("transcript_vars_objects_from_ensembl.py").read()) # sourced from create_transcript_objects_from_ensembl.py
except:                                                          # (dependant on downloaded files from ensembl website)
    print("could not load transcript_vars_objects_from_ensembl.py")

try:
    exec(open("transcript_vars_objects_from_sql.py").read()) # sourced from create_transcript_objects_from_sql.py
except:                                                      # (dependant on GTF and FASTA files which have been processed into a DB)
    print("could not load transcript_vars_objects_from_sql.py")


############
# FUNCTIONS:
############

#################################
def open_text_file(filename, rw):
#################################

    """
    This function opens a text file either for read or write.
    If any error is encountered, the program will abort (exit)
    
    Inputs: 2 strings, a filename and the open type ("r"ead/"w"rite)
    Output: File Handle
    """
    
    try:
        file_handle = open(filename, rw)
    except:
        print("Could not open",filename)
        exit()

    return file_handle


################
# Configuration:
################

"""
Set up as follows:

2 groups needed
 - as many transcripts as you like per group
 - as many mRNAs per transcript as you like (NB - If there are 3 junctions (i.e. 4 exons)
      in a transcript, asking for 1 mRNA will give you three RNA-Seq Reads - one per junction)
      

groupX = [
          [transcript1 (instance of Transcript class), #_mRNAs_to_generate (integer)],
          [transcript2 (instance of Transcript class), #_mRNAs_to_generate (integer)],
          .
          .
          .
          [transcriptN, #_mRNAs_to_generate]
         ]
"""


group1 = [
#          [PTF1A_ENST00000376504_4_001_002, 400000],
          [PHPT1_ENST00000371661_5_001_002_003_004_005, 400]
         ]

group2 = [
#          [PTF1A_ENST00000640579_1_001_002, 4250],
          [PHPT1_ENST00000247665_12_001_002_003, 1000]
         ]
         


groups = [group1, group2]

replicates_per_group = 5   # usually 5 (so five replicates per group)
read_length          = 50  # usually min 50 for RNA-Seq, between 28-32 for Ribo-Seq

# if we have read_lenth == 100 then left_offset = 50 and right_offset = 50 (for example)
# if we have read_lenth ==  99 then left_offset = 50 and right_offset = 49 (for example)
left_offset  = read_length // 2
right_offset = read_length // 2

if  left_offset + right_offset < read_length:
    left_offset += 1


############
# MAIN CODE:
############

"""
First we derive a list of junction offsets for each Transcript object/instance

If for example we have Exon Lengths of 22, 32, 40, 10 then the Junctions will be 22, 54, 94 (always one junction less than there are exons)

NB: A check is needed before populating the junction_list attribute of each Transcript object because the same Transcript object can appear more than once across the groups
    - without this check, every time the Transcript object is encountered, it will have more junctions added to its junction_list!
"""

current_offset = 0

for group in groups:
    for transcript, _ in group:   # _ is a "don't care" placeholder - we need it in order to skip over #_mRNAs_to_generate integer

        if  len(transcript.junction_list) == 0:
        
            current_offset = 0
            
            for each_exon_len in transcript.exon_len_list:
                current_offset += each_exon_len
                transcript.junction_list.append(current_offset)
            transcript.junction_list.pop()    # clip the last one (rightmost) off the list because we don't want the final NT in the transcript seq. to count as a junction


# sanity check to verify that Junctions are generated correctly from Exon Lengths
print("Exons Lengths -- Junction Offsets")
for group in groups:
    for transcript,total_reads in group:
        print(transcript.exon_len_list, " -- ", transcript.junction_list)
    print("\n")

##exit()


"""
Now we generate a .fastq file for each group
We'll end up with a total of 2 (groups) x replicates_per_group variable .fastq files
"""

curr_group = 0

for group in groups:

    curr_group += 1

    for curr_replicate in range(1, replicates_per_group + 1):

        filename_out = "synth_data_grp_" + str(curr_group) + "_replicate_" + str(curr_replicate) + ".fastq"
        output_file = open_text_file(filename_out, "w")
    
        for transcript, mRNAs_to_create in group:

            for i in range(1, mRNAs_to_create + 1):

                junction_count = 0

                for each_junction_offset in transcript.junction_list:
                    junction_count += 1
                    
                    start_pos = each_junction_offset - left_offset
                    end_pos   = each_junction_offset + right_offset

                    if  start_pos < 0:                             # have we gone too far left such that the RNA-Seq read starts before the transcript?
                        correction = abs(start_pos)
                        start_pos += correction
                        end_pos   += correction
                    if  end_pos > len(transcript.full_sequence):   # have we gone too far right such that the RNA-Seq read ends after the transcript?
                        correction = end_pos - len(transcript.full_sequence)
                        start_pos -= (correction + 1)
                        end_pos   -= (correction + 1)

                # Line 1
                    line = "@SYNTH_DATA " + "length=" + str(read_length) + " Read No: " + str(i) + " Junction: " + str(junction_count) + "/" + str(len(transcript.junction_list)) + " Xcript Name: " + transcript.name +" in Group: " + str(curr_group)
                    output_file.write(line)
                    output_file.write("\n")

                # Line 2
                    # we may look at "oscillating" the centre point of the read to invoke vague randomness (but read-fragments must include exon-exon junctions as it stands)
                    next_read = transcript.full_sequence[start_pos : end_pos]
                    output_file.write(next_read)
                    output_file.write("\n")

                # Line 3
                    line = "+SYNTH_DATA " + str(i) + " length=" + str(read_length) + " Phred Scores:"
                    output_file.write(line)
                    output_file.write("\n")

                # Line 4
                    line = "I" * read_length
                    output_file.write(line)
            ##        if  i != num_of_reads:
            ##            output_file.write("\n")  WHEN PUTTING GROUPS INTO LOOPS THIS CODE WILL NEED TO BE REWRITTEN TO CHECK FOR FINAL GROUP AND FINAL SAMPLE (NOT JUST FINAL SAMPLE)
                    output_file.write("\n")
                
        output_file.close()
