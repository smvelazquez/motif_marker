#!/Users/samanthavelazquez/miniconda3/bin/python

## Motif Marker v1.3:
## Samantha Velazquez
## 3/13/2020
## --------------------


import regex as re
import cairo
import math
import matplotlib
import argparse

def get_args():
     parser = argparse.ArgumentParser()
     parser.add_argument("-f", "--FASTA", help="the absolute FILEPATH for your FASTA file", required = True)
     parser.add_argument("-m", "--MOTIF", help="the absolute FILEPATH with your known motifs", required = True)
     return parser.parse_args()
args = get_args()

test_fasta = args.FASTA
test_motifs = args.MOTIF
motifs_list = []
orig_motifs = []

def switch_base(seq):
    """Takes a u/U and converts it to a regular expression in the form of [tT]"""
    invalid_bases = ['u', 'U']
    comp_dict = {'U': '[tT]', 'u':'[tT]'}
    bases=[]
    for base in seq:
        if base in invalid_bases:
            bases.append(comp_dict[base])
        else: bases.append(base)
    new_seq=bases
    return "".join(new_seq)

def ambiguous_bases(seq):
    """Function that creates regular expressions for every ambiguous bases.
        Y/y are now [tcTC]"""
    ambig_bases = ['y', 'Y']
    ambig_dict = {'Y':'[tcTC]', 'y': '[tcTC]'}
    bases=[]
    for base in seq:
        if base in ambig_bases:

            bases.append(ambig_dict[base])
        else: bases.append(base)

    return "".join(bases)

rgb_colors = [[0,0.8,0.8],[0,0,0.8],[0.6,0,0.6],[0.8,0,0],[0.25,0.25,0.25],[0.5,0.5,0.5],[0.6,0,0.3],[0.8,0.4,0],[0.8,0,0.8],[0,0.8,0]]
with open(test_motifs, "r") as fh:
    for line in fh:
        line = line.strip("\n")
        #print(line)
        orig_motifs.append(line)
        cleaned = switch_base(line) #change all of the u/U's to t/T's

        if re.finditer('Y', cleaned) or re.finditer('y', cleaned): ##search for any ambiguous bases
            new = ambiguous_bases(cleaned) #create the regular expression for them that fixes this
        else:
            continue
        #print(new, ':', cleaned)
        motifs_list.append(new)

width, height = 1000, 1000
#create the coordinates to display your graphic, desginate output
surface = cairo.SVGSurface("final_image.svg", width, height)

with open(test_fasta, "r") as file:
    records =0
    motifs = 0
    i=0
    original =0

    #writing out the legend
    context = cairo.Context(surface)
    context.move_to(10,400)
    context.set_source_rgb(0, 0, 0)
    context.show_text("Motif Key:")
    context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)


    for line in file:
        line = line.strip('\n')

        #create color dictionary
        c_index=0
        motif_color={}
        for mot in motifs_list:
            i+=1
            mot=mot.strip()
            motif_color[mot.strip()]=rgb_colors[c_index]
            x=20
            c_index+=1

        #pull out gene names for image
        if line.startswith('>'):
            header = line.split()
            gene = header[0]
            #print(gene)

            # write them next to their DNA strand
            context.set_font_size(10)
            context.set_source_rgb(0, 0, 0)
            context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            context.move_to(5,30*(records+1)+80)  # x,y coords; spacing graphs by 30px
            context.show_text(gene)

        ## remove any header lines
        if not line.startswith(">"):
            records +=1
            length = len(line) #store the length of the lines for use when drawing the picture

            #start drawing a line for every record, with the length of the sequence to scale
            context.move_to(50,30*(records+1)+50) # x,y coords; spacing graphs by 30px
            context.line_to(length+50,30*(records+1)+50) # Drawing line px length of sequence
            context.stroke()

            #find the exons
            for exon in re.finditer(r'[ATCG]+', line):
                exon_starts = exon.start()
                exon_end = exon.end()
                exon_coords = [exon_starts, exon_end]

                # Draw Rectangle to represent Exon
                context.rectangle(exon_coords[0]+50,(30*(records+1)+50)-10,exon_coords[1]-exon_coords[0],18)
                context.fill()

            #find the motifs and store them
            for motif in motifs_list:

                for match in re.finditer(motif, line, overlapped = True):

                    s = match.start()
                    e = match.end()
                    positions = [s,e]
                    motifs +=1

                    #set and color your motifs
                    context.set_source_rgb(motif_color[motif][0],motif_color[motif][1],motif_color[motif][2])
                    context.rectangle(positions[0]+50,(30*(records+1)+50)-10,positions[1]-positions[0],18)
                    context.fill()

    n=0
    for m in motif_color:
        ypos=len(orig_motifs)*100+10+n*20
        R = motif_color[m][0]
        G = motif_color[m][1]
        B = motif_color[m][2]
        context.set_source_rgb(R,G,B)
        context.rectangle(15,ypos,15,15)
        context.fill()
        context.stroke()
        context.move_to(40,ypos+10)
        motif_name=orig_motifs[n]
        context.show_text(motif_name)
        n+=1


surface.finish()
