import numpy as np
import pandas as pd

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from Bio.SubsMat import MatrixInfo as matlist



fasta_contents = open('/stanley/genetics/analysis/epi25/sc/wave5/misc/Homo_sapiens.GRCh38.pep.all.fa').read().strip().split('\n')
ensembl2sequence_dictionary = dict()
current_header = None
sequence_buffer = ''
for line in fasta_contents:
    if '>' in line:
        if current_header is not None:
            # current_header = current_header.split(".")[0] ### include version
            ensembl2sequence_dictionary[current_header] = sequence_buffer
        current_header = line.strip('>').split(' ')[0]
        sequence_buffer = ''
    else:
        sequence_buffer += line.strip()


fasta_contents = open('/stanley/genetics/analysis/epi25/sc/wave5/misc/uniprot_sprot.fasta').read().strip().split('\n')
uniprot2sequence_dictionary = dict()
current_header = None
sequence_buffer = ''
for line in fasta_contents:
    if '>' in line:
        if current_header is not None:
            current_header = current_header.split("|")[1]
            uniprot2sequence_dictionary[current_header] = sequence_buffer
        current_header = line.strip('>').split(' ')[0]
        sequence_buffer = ''
    else:
        sequence_buffer += line.strip()


def ensembl_to_sequence(ensembl):
    if ensembl in ensembl2sequence_dictionary:
        return ensembl2sequence_dictionary[ensembl]
    else:
        return str(np.nan)
    
    
def uniprot_to_sequence(uniprot):
    if uniprot in uniprot2sequence_dictionary:
        return uniprot2sequence_dictionary[uniprot]
    else:
        return str(np.nan)


# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def alignPident(align1, align2, useShorter=True):
    # Convert alignments to array, sum up identical AA
    identical = sum(np.array((list(align1)) == np.array(list(align2)))&(np.array(list(align1)) != "-"))
    # Figure out what to calculate percentage based off of
    # Shorter Sequence
    # Longer Sequence
    # Aligned Portion?
    # Total Alignment Length?
    if useShorter:
        totalLength = min([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
    else:
        totalLength = max([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
    # Return Percentage
    return float(identical)/totalLength
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def alignPositives(align1, align2, matrix=MatrixInfo.blosum62, useShorter=True):
    # Function for reading blosum score (since matrix uses non-symetrical (x, y) tuples as keys)
    def getBlosumScore(x, y):
        if (x, y) in matrix:
            val = matrix[(x, y)]
        elif (y, x) in matrix:
            val = matrix[(y, x)]
        else:
            val = 0
        return val
    # FUNCTION END
    # Calculate Blosum score over alilgnment
    blosumScores = [getBlosumScore(align1[i], align2[i]) for i in range(len(align1)) if not (align1[i] == "-" or align2[i] == "-")]
    # Sum Positives
    positive = sum(np.array(blosumScores) > 0)
    # Figure out what to calculate percentage based off of
    # Shorter Sequence
    # Longer Sequence
    # Aligned Portion?
    # Total Alignment Length?
    if useShorter:
        totalLength = min([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
    else:
        totalLength = max([sum(np.array(list(align1)) != "-"), sum(np.array(list(align2)) != "-")])
    # Return Percentage
    return float(positive)/totalLength
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def alignCoverage(align1, align2, matrix=MatrixInfo.blosum62, useShorter=True):
    # Sum up non-gapped regions of alignment
    non_gaps = len(align1) - sum([align1[i] == "-" or align2[i] == "-" for i in range(len(align1))])
    # Divide by sequence length of align1
    return float(non_gaps)/len(align1.replace("-", ""))
# FUNCTION END

# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def getBlosumScore(x, y, matrix=MatrixInfo.blosum62):
    if (x, y) in matrix:
        val = matrix[(x, y)]
    elif (y, x) in matrix:
        val = matrix[(y, x)]
    else:
        val = 0
    return val
# FUNCTION END


def matrix_get(aa1, aa2):
    matrix = matlist.blosum62
    if (aa1, aa2) in matrix:
        val = matrix[(aa1, aa2)]
    elif (aa2, aa1) in matrix:
        val = matrix[(aa2, aa1)]
    else:
        val = 0
    return val


# Alignment Helper Functions (Copied from Rice Eclair Predictions Notebook July 17, 2018)
def NWSeqAlignment(s1, s2, matrix=MatrixInfo.blosum62, useShorter=True, best=False):
    keep = None
    for a in pairwise2.align.globalcs(s1, s2, matrix_get, -10, -0.5):
        a1, a2, score, begin, end = a
        r = dict([("Align1", a1), ("Align2", a2), ("Score", score), ("Begin", begin), ("End", end), ("Pident", alignPident(a1, a2, useShorter=useShorter)), ("Positives", alignPositives(a1, a2, matrix=matrix, useShorter=useShorter)), ("Coverage1", alignCoverage(a1, a2, matrix=matrix, useShorter=useShorter)), ("Coverage2", alignCoverage(a2, a1, matrix=matrix, useShorter=useShorter))])
        align_string = pd.DataFrame([list(r["Align1"]), list(r["Align2"])]).T
        def do(z):
            x,y = z
            if(x == y):
                return "|"
            if(x == "-" or y == "-"):
                return " "
            else:
                return ["-", "+"][getBlosumScore(x, y, matrix=matrix) > 0]
        # FUNCTION END
        align_string[2] = list(zip(align_string[0], align_string[1]))
        align_string["C"] = align_string[2].apply(do)
        align_string = "".join(align_string["C"].values)
        r["Alignment"] = align_string
        
        if(keep == None or r["Pident"] >= keep["Pident"]):
            keep = r
        if(not best):
            break
    return keep
# FUNCTION END

# Save mapping between MSU Seq positions and Uniprot Seq positions
def map_positions(align1, align2):
    pos1 = 1
    pos2 = 1
    rows = []
    d1 = {}
    d2 = {}
    for i in range(len(align1)):
        AA1 = align1[i]
        AA2 = align2[i]
        if(align1[i] != "-" and align2[i] != "-"):
            d1[pos1] = pos2
            d2[pos2] = pos1
            tmppos1 = pos1
            tmppos2 = pos2
            pos1 += 1
            pos2 += 1
        else:
            if(align1[i] != "-"):
                d1[pos1] = None
                tmppos1 = pos1
                tmppos2 = ""
                pos1 += 1
            if(align2[i] != "-"):
                d2[pos2] = None
                tmppos1 = ""
                tmppos2 = pos2
                pos2 += 1
        rows.append((tmppos1, tmppos2, AA1, AA2))
    return d1, d2, rows
# FUNCTION END
