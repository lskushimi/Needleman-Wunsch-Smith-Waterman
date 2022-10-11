#!/usr/bin/env python

import sys
import argparse
from turtle import distance
import numpy as np

base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3
ALIGN_GLOBAL, ALIGN_SEMIGLOBAL, ALIGN_LOCAL = 0, 1, 2

def seqalignDP(seq1, seq2, subst_matrix, gap_penalty, align_type):
    """
    Return the score of the optimal Needdleman-Wunsch alignment for seq1
    and seq2.
    Note: gap_penalty should be positive (it is subtracted)
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]


    # YOUR CODE HERE
    # Fill in the dynamic programming tables F and TB
    # Hints: The first row and first column of the table F[i][0] and F[0][j]
    # contain dummy values.
    #  (see for illustration Durbin p.21, Figure 2.5, but be careful what you
    #   think of as rows and what you think of as columns)
    #  Hence, the bases corresponding to F[i][j] are actually seq1[i-1] and
    #  seq2[j-1].
    #  Use the dictionary base_idx to convert from the character to an index to
    #   look up entries of the substitution matrix.
    #  To get started, you can complete and run the algorithm filling in only
    #    F, and then figure out how to do TB.
    # YOU NEED TO handle three different values for align_type: ALIGN_GLOBAL, ALIGN_SEMIGLOBAL, ALIGN_LOCAL
    # (defined above)

    matchTemp = 0
    gapTempX = 0
    gapTempY = 0
    x = len(seq1)
    y = len(seq2)

    #initialization of global alignment
    if align_type == ALIGN_GLOBAL:
        for i in range(len(seq1)+1):
            F[i][0] = -gap_penalty * i

        for j in range(len(seq2)+1):
            F[0][j] = -gap_penalty * j

    #initialization of semiglobal
    elif align_type == ALIGN_SEMIGLOBAL:
        if len(seq1) >= len(seq2):
            for i in range(len(seq1) + 1):
                F[i][0] = 0
            for j in range(len(seq2)+1):
                F[0][j] = -gap_penalty * j
        else:
            for j in range(len(seq2) + 1):
                F[0][j] = 0
            for i in range(len(seq1)+1):
                F[i][0] = -gap_penalty * i
    
    #initialization of local
    elif align_type == ALIGN_LOCAL:
        for i in range(len(seq1) + 1):
            F[i][0] = 0
        for j in range(len(seq2) + 1):
            F[0][j] = 0


    #Recursion is the same for global and semi-global
    if align_type == ALIGN_GLOBAL or align_type == ALIGN_SEMIGLOBAL:
        for i in range(1,len(seq1)+1):
            TB[i][0] = PTR_GAP2
            for j in range(1,len(seq2)+1):
                if seq1[i-1] == seq2[j-1]:
                    matchTemp = F[i-1][j-1] + subst_matrix[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]
                else:
                    matchTemp = F[i-1][j-1] + subst_matrix[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]
                gapTempY = F[i-1][j] - gap_penalty
                gapTempX = F[i][j-1] - gap_penalty
                F[i][j] = max(matchTemp,gapTempX,gapTempY)
                if F[i][j] == matchTemp:
                    TB[i][j] = PTR_BASE
                elif F[i][j] == gapTempX:
                    TB[i][j] = PTR_GAP1
                elif F[i][j] == gapTempY:
                    TB[i][j] = PTR_GAP2
                TB[0][j] = PTR_GAP1

    #different recursion for local alignment
    elif align_type == ALIGN_LOCAL:
        for i in range(1,len(seq1)+1):
            TB[i][0] = PTR_GAP2
            for j in range(1,len(seq2)+1):
                if seq1[i-1] == seq2[j-1]:
                    matchTemp = F[i-1][j-1] + subst_matrix[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]
                else:
                    matchTemp = F[i-1][j-1] + subst_matrix[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]
                gapTempY = F[i-1][j] - gap_penalty
                gapTempX = F[i][j-1] - gap_penalty
                F[i][j] = max(matchTemp,gapTempX,gapTempY,0)
                if F[i][j] == matchTemp:
                    TB[i][j] = PTR_BASE
                elif F[i][j] == gapTempX:
                    TB[i][j] = PTR_GAP1
                elif F[i][j] == gapTempY:
                    TB[i][j] = PTR_GAP2
                elif F[i][j] == 0:
                    TB[i][j] = PTR_NONE
                TB[0][j] = PTR_GAP1

    #0 out after finding max value in last row or column for semi-global
    if align_type == ALIGN_SEMIGLOBAL:
        if len(seq1) > len(seq2):
            Fprime = np.array(F)
            maxIdx = np.where(Fprime[:,len(seq2)] == max(Fprime[:,len(seq2)]))
            for i in range(maxIdx[0][0] + 1, len(seq1)+1):
                for j in range(len(seq2)+1):
                    F[i][j] = 0
                    TB[i][j] = 0
        #change which index used for scoring
        while F[x][y] == 0:
            if F[x-1][y] == 0:
                y = y-1
            elif F[x][y-1] == 0:
                x = x-1

    if align_type == ALIGN_LOCAL:
        Fprime = np.array(F)
        maxVal = np.max(Fprime)
        maxIdx = np.where(Fprime == maxVal)
        x = maxIdx[0][0]
        y = maxIdx[1][0]

        #zero out values not relevant to the max
        for i in range(maxIdx[0][0] + 1, len(seq1)+1):
            for j in range(len(seq2) + 1):
                F[i][j] = 0
                TB[i][j] = 0

        for i in range(1,len(seq1) + 1):
            for j in range(maxIdx[1][0] + 1, len(seq2)+1):
                    F[i][j] = 0
                    TB[i][j] = 0
        for i in range(len(seq1)):
            TB[i][0] = 0
        for j in range(len(seq2)):
            TB[0][j] = 0

    print(np.matrix(F))
    print(np.matrix(TB))

    return F[x][y], F, TB


def traceback(seq1, seq2, TB):

    # ONLY VALID FOR GLOBAL ALIGNMENTS
    # You will need to modify this for semiglobal and local

    s1 = ""
    s2 = ""

    i = len(seq1)
    j = len(seq2)

    TBprime = np.array(TB)

    #for semi-global
    if TB[i][j] == 0 and (TB[i-1][j] != TB[i][j-1] == 0):

        while TB[i][j] == 0:
            if TB[i-1][j] == 0:
                s1 = '-' + s1
                s2 = seq2[j-1] + s2
                j = j-1

            elif TB[i][j-1] == 0:
                s1 = seq1[i-1] + s1
                s2 = '-' + s2
                i = i-1

    #for local
    if TB[i][j] == 0 and TB[i-1][j] == 0 and TB[i][j-1] == 0:
        startPoint = np.where(TBprime != 0)
        startPoint = np.flip(startPoint)
        max_i = startPoint[1][0]
        max_j = startPoint[0][0]
        i = max_i
        j = max_j


    # while TB[i][j] == 0:
    #     for x in reversed(range(1,i+1)):
    #         for y in reversed(range(1,j+1)):
    #             if TB[x][y] != 0:
    #                 i = x
    #                 j = y

        #print(x,y)
        #print(TB[x][y])


    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i = i - 1
            j = j - 1
        elif TB[i][j] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j = j - 1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i = i - 1
        else:
            assert False

    distanceScore = 0

    S = [
    # A   G   C   T
    [ 3, -1, -2, -2],  # A
    [-1,  3, -2, -2],  # G
    [-2, -2,  3, -1],  # C
    [-2, -2, -1,  3]   # T
    ]    
    
    for i in range(len(s1)-1):
        if s1[i] == s2[i]:
            distanceScore = distanceScore + 0
        elif s1[i] == '-' and s2[i] != '-':
            distanceScore = distanceScore + 4
        elif s1[i] != '-' and s2[i] == '-':
            distanceScore = distanceScore + 4
        elif s1[i] != '-' and s2[i] != '-' and s1[i] != s2[i]:
            distanceScore = distanceScore - S[base_idx[s1[i]]][base_idx[s2[i]]]

    print('Distance score is: ' + str(distanceScore))

    return s1, s2


def readSeq(filename):
    """Reads in a FASTA sequence. Assumes one sequence in the file"""
    seq = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.rstrip().upper())

    return "".join(seq)

# Substituation matrix and gap_penalty
S = [
    # A   G   C   T
    [ 3, -1, -2, -2],  # A
    [-1,  3, -2, -2],  # G
    [-2, -2,  3, -1],  # C
    [-2, -2, -1,  3]   # T
]
gap_penalty = 4


def main():
    # parse command line
    parser = argparse.ArgumentParser(usage='%(prog)s --type [alignment_type] <fasta1> <fasta2>')
    
    
    parser.add_argument('--type', default=ALIGN_GLOBAL, 
                        choices=[0,1,2],
                        type=int, 
                        help="the type of alignment: 0=global (default), 1=semiglobal, 2=local")

    parser.add_argument('fasta1', help="fasta file 1")
    parser.add_argument('fasta2', help="fasta file 1")
    
    args=parser.parse_args()
    
    
    align_type = args.type    
    file1 = args.fasta1
    file2 = args.fasta2

    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    score, F, TB = seqalignDP(seq1, seq2, S, gap_penalty, align_type)

    print("Score: {0}".format(score))

    s1, s2 = traceback(seq1, seq2, TB)
    print(s1)
    print(s2)

if __name__ == "__main__":
    main()
