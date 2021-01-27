import sqlite3
import collections # for deque

def connect_to_sqlite_db(db_name):

    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    return c

#####################
##  CONFIG VARIABLES:
#####################

##file_name = "gencode.v26.annotation.gtf"
##db_name   = "fasta_gtf_db.db"

file_name = "gencode.v35.primary_assembly.annotation.gtf"
db_name   = "fasta_gtf_v35.db"


############
# MAIN CODE:
############

try:
    infile = open(file_name, "r")
except:
    print("input file error on open...halting process")
    exit()


line_count = 0

conn = sqlite3.connect(db_name)
c = conn.cursor()

for line in infile:

    if  line == "" or line[0:3] != "chr":   # do not change the order of the IF condition! Left is evaluated before right is looked at in Python!
        continue

    split_string = line[:-1].split("\t")    # split the line based on tab, and also remove last char as it will be CR/LF

    if  split_string[2] == "exon":

        line_count += 1

        extra_info = split_string[8].split(";")  # split the extra info (column 9) into component parts

        CoordFrom  = split_string[3]
        CoordTo    = split_string[4]
        Chrom      = split_string[0][3:] # assuming all chromosome columns start with chr and are in format chr1, chr19, chrM etc
        Strand     = split_string[6]
        
        for each_name_value_pair in extra_info:
            each_name_value_pair = each_name_value_pair.strip()
            
##            print(each_name_value_pair)
##            print()

            split_NV_pair = each_name_value_pair.split(" ")

##            print(split_NV_pair)

            if  split_NV_pair[0] == "gene_name":
                GeneName = split_NV_pair[1].replace('"', '')
            elif split_NV_pair[0] == "transcript_id":
                XcriptID = split_NV_pair[1].replace('"', '')
            elif split_NV_pair[0] == "exon_number":
                ExonNo = split_NV_pair[1].replace('"', '')
                ExonNo = ("000" + ExonNo)[-3:]  # we need ExonNo to be consecutive string, so left pad it with zeros to make 3-len string. Max ExonNo==999
            elif split_NV_pair[0] == "transcript_type":
                XcriptType = split_NV_pair[1].replace('"', '')
            elif split_NV_pair[0] == "gene_id":
                GeneID = split_NV_pair[1].replace('"', '')
            elif split_NV_pair[0] == "exon_id":
                ExonID =split_NV_pair[1].replace('"', '')
            elif split_NV_pair[0] == "transcript_name":
                XcriptName = split_NV_pair[1].replace('"', '')

        row = (GeneName, XcriptID, ExonNo, CoordFrom, CoordTo, XcriptType, GeneID, ExonID, XcriptName, Chrom, Strand)
##        print(row)
        
        c.execute('insert into EXONS_GTF values (?,?,?,?,?,?,?,?,?,?,?)', row) # row is a tuple of elements. One element for each "?" placeholder

        if  line_count >= 10000:   # we flush out the buffer after this many INSERTs by doing a COMMIT of the data
            conn.commit()
            line_count = 0    # reset this so we're ready to package up another 100 rows to commit
        
conn.commit() # commit for the last set of rows (when less than 10000 left - should be ok if nothing to commit too, if file row-count // 100 == 0)
conn.close()
