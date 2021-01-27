import sqlite3
import collections # for deque


###########
# DB SETUP:
###########

#db_name = 'fasta_gtf_v26.db'
db_name = 'fasta_gtf_v35.db'

conn = sqlite3.connect(db_name)
c = conn.cursor()

############
# MAIN CODE:
############

#file_name = "GRCh38.p10.genome.fa"   # hg38 (should all be the same - fasta only changes when hgNN changes, not when Version Number changes)
file_name = "GRCh38.v35.primary_assembly.genome.fa"

try:
    infile = open(file_name, "r")
except:
    print("input file error on open...halting process")
    exit()


end_reached = False
deq_all_lines = collections.deque()

commit_line_count = 0
first_chrom = True

for line in infile:

    if  line == "":
        continue            # we don't care about blank lines - go back to the "for" loop for the next iteration

    if  line[0] == ">":     # we've moved on to a "new" (or the first) chromosome

        if  line[0:4] != ">chr":
            end_reached = True

        if  first_chrom:

            split_header = line[:-1].split(" ")    # the header format is ">chr1 1" so "split it on space"... 
            chrom_no = split_header[1]             # ...and take the part after the space (i.e. Chromosome number: 1, 2, 14, X, Y, M etc)

            first_chrom = False
            
            continue        # no nucleotides to process when we encounter the first chromosome header - so go back to the "for" loop for the next iteration
        
        else: # we've changed chromosome number, and we're not dealing with the first chromosome, so we need to write full chromosome (stored in deq_all_lines) to DB
         
            fasta_1000_line = collections.deque()    # initialise a new deque (this will take 1000 chars at a time from the full-chromosome deque deq_all_lines)
            row_count = 1                            # restart row_count for each chromosome 

            while len(deq_all_lines) > 0:

                fasta_1000_line.extend(deq_all_lines.popleft()) # deque.extend: takes an iterable (eg list/tuple/string)...
                                                                # ...and adds each element of the iterable to the list one at a time
                                                                
                if  len(fasta_1000_line) > 999 \
                or  len(deq_all_lines) == 0:                    # if we have moved 1000 nucleotides to fasta_1000_line (ready to write to DB)...# ...or we have no more nucleotides left in the full-chromosome deq_all_lines
                    
                    str_fasta_1000_line = ''.join(fasta_1000_line) # ...this line might be s.l.o.w.! (convert DB-relevant deque to string)
                    
                    rec_fasta_1000_line = (chrom_no, str(row_count), str_fasta_1000_line,) # create a tuple, ready for DB-write

                    commit_line_count += 1
                    c.execute('insert into FASTA_1000 values (?,?,?)', rec_fasta_1000_line) # rec_fasta_1000_line is a tuple of elements. One element for each "?" placeholder
                    if  commit_line_count > 10000:
                        conn.commit()
                        commit_line_count = 0
                    
                    row_count += 1
                    fasta_1000_line = collections.deque()       # blank out the deque which does the 1000-nucleotide writes

            deq_all_lines = collections.deque() # reset deq_all_lines to be an empty deque
            
        split_header = line[:-1].split(" ")  # the header format is ">chr1 1" so "split it on space"... 
        chrom_no = split_header[1]           # ...and take the part after the space (i.e. Chromosome number: 1, 2, 14, X, Y, M etc)

    else:
        deq_all_lines.extend(line[:-1]) # add the input-file"line" to the end of the deque (don't include last char of "line" as it is CR/LF)

    if  end_reached: break

conn.commit() # commit for the last set of rows (if anything left uncommitted - especially at end of final chromosome!)
conn.close()  # close the connection and we're all done!
