import matplotlib.pyplot as plt
import numpy as np
import copy   # for deepcopy of lists

"""

======================
Usage Notes:

- This program relies on an input file called: input_filename = "bar_chart_config.txt"
- ALL TAKEN FROM A GENCODE GTF (version V35 of HG38 at time of writing)!!
- Every transcript must be followed by a # symbol
- Fields are tab-separated
- Any transcript can be used, including parts of transcripts.
- Any exons can be included/omitted, but the included exons should be ascending
- Each entry represents: GeneName - TranscriptId - ExonNo - CoordFrom - CoordTo

An example:

PHPT1	ENST00000247665.11	001	136849094	136849590
PHPT1	ENST00000247665.11	002	136850013	136850137
PHPT1	ENST00000247665.11	003	136850755	136851022
#
PHPT1	ENST00000371661.5	001	136848724	136848780
PHPT1	ENST00000371661.5	002	136849111	136849590
PHPT1	ENST00000371661.5	003	136850013	136850137
PHPT1	ENST00000371661.5	004	136850506	136850560
PHPT1	ENST00000371661.5	005	136850755	136851022
#
PHPT1	ENST00000462205.5	001	136849425	136849590
PHPT1	ENST00000462205.5	002	136850013	136850137
PHPT1	ENST00000462205.5	003	136850503	136850560
PHPT1	ENST00000462205.5	004	136850755	136851024
#
PHPT1	ENST00000463215.1	001	136849779	136850137
PHPT1	ENST00000463215.1	002	136850755	136851022
#
PHPT1	ENST00000492540.5	001	136849091	136849590
PHPT1	ENST00000492540.5	002	136850013	136850560
PHPT1	ENST00000492540.5	003	136850755	136851022
#
PHPT1	ENST00000497413.1	001	136849426	136849590
PHPT1	ENST00000497413.1	002	136850013	136850137
PHPT1	ENST00000497413.1	003	136850506	136850560
PHPT1	ENST00000497413.1	004	136850755	136851027
#

- The Horizontal Bar Chart produced will allow easy visualisation of different transcripts in the same gene,
  as we are interested in how the exons and introns overlap.

- The transcript with the lowest Coord-From value in the transcripts present in "bar_chart_config.txt"
  will be given an offset of 0 and be the leftmost transcript in the horizontal bar chart.

- A summary of the offsets of the exons of each transcript is also output by the program, along with a
  summary of what has been included in "bar_chart_config.txt". All useful for experiment documentation.

=======================

Programming Flow/Technique notes:

- fill this in later, following step-by-step how the program works (also see code comments)

=======================

Exon To and Froms (no introns yet) - here we have a list, containing 2 sublists,
each contain lists of pairs of To-From coords
[[[136850013, 136850137],
  [136850506, 136850560],
  [136850755, 136851022]],
 [[136849111, 136849590],
  [136850013, 136850137],
  [136850506, 136850560],
  [136850755, 136851022]]]

136849111
[[[136850013, 136850137], [136850506, 136850560], [136850755, 136851022]], [[136849111, 136849590], [136850013, 136850137], [136850506, 136850560], [136850755, 136851022]]]
[[[902, 1026], [1395, 1449], [1644, 1911]], [[0, 479], [902, 1026], [1395, 1449], [1644, 1911]]]
"""


# The Xcript class will hold xcript instances which will be printed out at the end of the run (for documentation)
# This class is not used in any of the processing to create the matplotlib bar chart - it's only used for the summary print
class Xcript:
    def __init__(self, name, gene, xcript_id):
        self.name          = name
        self.gene          = gene
        self.xcript_id     = xcript_id
        self.exon_list     = []

list_of_xcript_objects = []

input_filename = "bar_chart_config.txt"
             
try:
    infile = open(input_filename, "r")
except:
    print("input file error on open...halting process")
    exit()


list_of_xcripts = []

curr_xcript = 0
list_of_xcripts.append([])

for line in infile:
    if  line[0] != "#":

        split_string = line[:-1].split("\t")

        list_of_xcripts[curr_xcript].append([int(split_string[3]), int(split_string[4])] )

############## BEGIN: Xcript Class / Object Processing ###############
        object_found = False                        
        for each_xcript in list_of_xcript_objects:   # scan through the existing Xcript objects...
            if  each_xcript.name == curr_xcript + 1: # ...if we find one with a ".name" match the current xcript count...
                object_found = True                  # ...stop executing the loop so we can add the new exon to the correct object
                break
        if  (not object_found):                      # create a new Xcript object if we don't have one for this curr_xcript yet
            new_xcript = Xcript(curr_xcript + 1, split_string[0], split_string[1]) # Transcript No (+1 as we don't want 0), Gene, Xcript_Id
            new_xcript.exon_list.append(split_string[2])                           # add new exon to new Xcript object
            list_of_xcript_objects.append(new_xcript)                              # append the new object to the list_of_xcript_objects
        else:
            each_xcript.exon_list.append(split_string[2])                          # add new exon to existing Xcript object
############## END: Xcript Class / Object Processing ###############            
    else:
        list_of_xcripts.append([])
        curr_xcript += 1

if  line[0] == "#":       # We finished on a hash meaning we have an empty list...
    list_of_xcripts.pop() # ...remove it!

lowest_start_coord = min([first_elem[0][0] for first_elem in list_of_xcripts]) # use list comprehension to get the lowest "From" Coord - this offset will be 0

###########################
## insert intron to-from coords (initialise their FROM-TO values to [-1, -1]:

exon_intron_coord_from_to_list = []

for i in range(0, len(list_of_xcripts)):
    
    exon_intron_coord_from_to_list.append([])
    
    for j in range(0, len(list_of_xcripts[i])):

        exon_intron_coord_from_to_list[i].append(list_of_xcripts[i][j])
        exon_intron_coord_from_to_list[i].append([-1,-1])

for i in range(0, len(exon_intron_coord_from_to_list)): # remove all the "last" [-1,-1] entries because we don't want our transcripts to finish on an intron
    exon_intron_coord_from_to_list[i].pop()


###########################
## overwrite all the [-1,-1] entries (i.e. the introns!) with ["coord of previous exon to + 1", "coord of next exon to - 1"]

for i in range(0, len(exon_intron_coord_from_to_list)):
    for j in range(0, len(exon_intron_coord_from_to_list[i])):
        if  exon_intron_coord_from_to_list[i][j] == [-1,-1]:
##            print("intron found")
            exon_intron_coord_from_to_list[i][j][0] = exon_intron_coord_from_to_list[i][j-1][1] + 1
            exon_intron_coord_from_to_list[i][j][1] = exon_intron_coord_from_to_list[i][j+1][0] - 1


###########################
## create a new list with the same exon & intron coords, but subtract the "lowest" exon-from-coord from ALL values to create a base-zero list of values
#   i.e. we'll "transform" it to be base-zero
#   (NB: deep copy needed obviously - how dreadfully poor is Python in this area?!)
lowest_start_coord = min([first_elem[0][0] for first_elem in exon_intron_coord_from_to_list])



exon_intron_coord_from_to_list_base_zero = copy.deepcopy(exon_intron_coord_from_to_list) # deep copy needed obviously (how dreadfully poor is Python in this area?!)

for i in range(0, len(exon_intron_coord_from_to_list_base_zero)):
    for j in range(0, len(exon_intron_coord_from_to_list_base_zero[i])):
        for k in range(0,2):
            exon_intron_coord_from_to_list_base_zero[i][j][k] -= lowest_start_coord


###########################
## ALL transcripts, including the one which now starts at "coord zero", need to have a "padding" area added. This is neither an intron nor an exon, and
##  it's needed for matplotlib to allow horizontal bars to "float" to the right of the (LHS of the) zero-based exon
##
## The transcript with the zero-based exon will have [0,0] inserted (to position 0) as its padding
## All other transcripts will have [0, n] inserted (to position 0) as its padding, (where n is first-exon-from MINUS 1)

for i in range(0, len(exon_intron_coord_from_to_list_base_zero)):
    if  exon_intron_coord_from_to_list_base_zero[i][0][0] != 0:
        n = exon_intron_coord_from_to_list_base_zero[i][0][0] - 1
        exon_intron_coord_from_to_list_base_zero[i].insert(0, [0, n])
    else:
        exon_intron_coord_from_to_list_base_zero[i].insert(0, [0,0])

############################
## Because matplotlib is fussy with the stacked horizontal bars (hbar), we need to have the same number of elements in each transcript
## We can define an "element" as the padding, exons, or introns.
## So we need to find the longest transcript, and APPEND the shorter ones with "empty" elements to pad them out to the same length as the longest transcript

# find length of longest xcript
all_lengths = []

for i in exon_intron_coord_from_to_list_base_zero:
    all_lengths.append(len(i))

longest_xcript_len = max(all_lengths)


# right-pad all xcripts that are shorter than the max to be the same length as the max (with [n,n] so there's no visual effect in matplotlib, where n = TO-COORD of last exon)

for i in exon_intron_coord_from_to_list_base_zero:

    if  len(i) < longest_xcript_len:
        last_exon_to = i[-1][1]    # we want the right-paired element in the last sub-list of this transcript

        for j in range(0, (longest_xcript_len - len(i))): # add as many extra [n,n] elements as we need to bring this transcript up to the same length as the longest transcript
            i.append([last_exon_to,last_exon_to])
            
################################
## Next we want to create a list for each transcript which contains the width of each element (where an element is padding, exon or intron)

exon_intron_widths_list = []
for i in range(0,len(exon_intron_coord_from_to_list_base_zero)):
    exon_intron_widths_list.append([])

for i in range(0, len(exon_intron_coord_from_to_list_base_zero)):
    for from_to_coords in exon_intron_coord_from_to_list_base_zero[i]:
        exon_intron_widths_list[i].append(from_to_coords[1] - from_to_coords[0])


#################################
## The widths of each element of each transcript is now stored in exon_intron_widths_list
##
## The next thing we need to do is rearrange the data from lists of transcripts to lists of offset-element lengths

elements_and_widths = []
for i in range(0, len(exon_intron_widths_list[0])): # any of the transcript lengths will do since they're all the same
    elements_and_widths.append([])


for i in range(0, len(exon_intron_widths_list[0])):  # any of the transcript lengths will do since they're all the same
    for j in range(0, len(exon_intron_widths_list)): # "how many transcripts" length
        elements_and_widths[i].append(exon_intron_widths_list[j][i])


#################################
## We can therefore draw the "relative" horizontal stacked bars

xcript_no   = []
offset_list = []

for i in range(1, len(elements_and_widths[0]) + 1):
    xcript_no.append(i)
    offset_list.append(0)


if  len(exon_intron_coord_from_to_list_base_zero) < 10:
    calc_height = len(exon_intron_coord_from_to_list_base_zero)/10
else:
    calc_height = 0.8
print(calc_height)
calc_height = 0.9

color_exon_intron = ["#ff0000", "#00ff00", "#ffffff"]  # red (introns), green (exons), white (left-padding)

# let's set up the padding first:
current_color = 2
plt.barh(xcript_no, elements_and_widths[0], height = calc_height, left = offset_list, color=color_exon_intron[current_color])


current_color = 1
offset_list = list(np.add(offset_list, elements_and_widths[0]))

for i in range(1, len(elements_and_widths)): # start at 1 because padding is at pos 0 and has already been done

    plt.barh(xcript_no, elements_and_widths[i], height = calc_height, left = offset_list, color=color_exon_intron[current_color])
    offset_list = list(np.add(offset_list, elements_and_widths[i]))
    current_color = abs(current_color - 1) # set the color to the opposite of whatever it was (oscillate between 0 and 1, i.e. red and green, i.e. intron and exon)
    
# each transcript number will be represented on the y-axis. We don't want decimal values to appear on the scale shown on the y-axis, so we
# can tell matplotlib to only show the transcript numbers (ranging from 1 to however many transcripts we have) by giving it a list of integers

transcript_numbers = range(1, len(exon_intron_coord_from_to_list_base_zero) + 1) # produces a list, e.g. [1, 2, 3, 4, 5, 6]

plt.yticks(transcript_numbers)
plt.ylabel("Transcript Number")
plt.xlabel("Nucleotide Offset - from coord:" + str(lowest_start_coord))
plt.legend(labels=[None, 'Exon', 'Intron']) # None for whitespace padding, then Exon and Intron as we added them in that order
plt.title("Intron / Exon overlaps")

plt.legend(labels=[None, 'Exon', 'Intron'], bbox_to_anchor=(0.80,0.99)) # None for whitespace padding, then Exon and Intron as we added them in that order
plt.tight_layout()
###################################
# OUTPUT REPORT:
###################################

###################################
# Show Transcripts (from Coord 0 to make it easier to read!), both Introns and Exons
#  PLUS whitespace and end-padding (to-from coords are identical):

print("##################################")
print("##################################")
print("Output Report:")
print("")
print("Nucleotide Offset - From Coordinate:", lowest_start_coord)
print("")

for i in range(0, len(exon_intron_coord_from_to_list_base_zero)):
    print("Transcript: ", i+1)
    for j in range(0, len(exon_intron_coord_from_to_list_base_zero[i])):
        print(exon_intron_coord_from_to_list_base_zero[i][j])

print("")
print("")

###################################
### PRINT FORMATTED TRANSCRIPTS - 3 COLUMNS PER ROW, EACH COL IS A TRANSCRIPT:
###  Whitespace and end-padding not shown

total_xcripts = len(exon_intron_coord_from_to_list_base_zero)
curr_xcript = 0
exon_string   = "Exon:    "
intron_string = "Intron:  " 
curr_component = exon_string

start_at_xcript = 0
stop_at_xcript  = 0

while stop_at_xcript < total_xcripts:
    
    stop_at_xcript  += 3
    if  stop_at_xcript > total_xcripts:
        stop_at_xcript = total_xcripts

    for i in range(start_at_xcript, stop_at_xcript):
        out_string = "Transcript: " + str(i+1) + " "*20
        out_string = out_string[0:27]
        print(out_string, end = "")
    print("")

    for j in range(1, len(exon_intron_coord_from_to_list_base_zero[0])):   # range from 1 as we want to skip the whitespace (any xcript len will do as all are same)
        for i in range(start_at_xcript, stop_at_xcript):
            if  exon_intron_coord_from_to_list_base_zero[i][j][0] != exon_intron_coord_from_to_list_base_zero[i][j][1]:
                out_string = curr_component + str(exon_intron_coord_from_to_list_base_zero[i][j]) + " "*20
                out_string = out_string[0:27]
                print(out_string, end = "")
            else:
                print(" "*27, end = "")
        if  curr_component == exon_string:
            curr_component = intron_string
        else:
            curr_component = exon_string
        print("")

    print("")
    curr_component = exon_string
    start_at_xcript += 3

## now we print out the summary of the config file:
print("Transcript Number ", "Gene Name ", "Transcript Id        ", "Exons")
for i in list_of_xcript_objects:
    out_string = " "*3 + str(i.name) + " "*17
    out_string = out_string[0:17]
    print(out_string, end = "")
    
    out_string = " "*2 + i.gene + " "*10
    out_string = out_string[0:12]
    print(out_string, end = "")

    out_string = " "*1 + i.xcript_id + " "*10
    out_string = out_string[0:22]
    print(out_string, end = "")

    out_string = " "*1 + str(i.exon_list) + " "*10
    out_string = out_string
    print(out_string, end = "")
    print()
        
##################################

plt.show() # show the plot at the very end, so that the Output Report can be generated first and both can be seen together
