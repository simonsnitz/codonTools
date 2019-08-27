from Bio import pairwise2 #allows us to do pairwise alignment on our proteins
from Bio.pairwise2 import format_alignment

from Bio.SubsMat import MatrixInfo as matlist #gives us access to blosum62 for alignment

from Bio.Seq import Seq #allows us to go from nt sequence to aa sequence
from Bio.Alphabet import generic_dna

def codon_optimize(rawinseq1, rawinseq2):


    rawinseq1nt = rawinseq1.upper().rstrip()
    for char in rawinseq1nt:
        if char != 'A' and char != 'T' and char != 'C' and char != 'G':
            print('\n'+'Error: input contains non-DNA bases(A,T,C,G)'+'\n'+'Continue running? Y or N')
            choice = raw_input()
            if choice != 'Y' and choice != 'y':
                quit()
    inseq1nt=Seq(rawinseq1nt, generic_dna)



    rawinseq2nt = rawinseq2.upper().rstrip()
    for char in rawinseq2nt:
        if char != 'A' and char != 'T' and char != 'C' and char != 'G':
            print ('\n'+'Error: input contains non-DNA bases(A,T,C,G)'+'\n'+'Continue running? Y or N')
            choice = raw_input()
            if choice != 'Y' and choice != 'y':
                quit()
    inseq2nt=Seq(rawinseq2nt, generic_dna)

    inseq1aa=inseq1nt.translate(to_stop=True) #note: this will translate until the stop codon so if you've got a bunch of stop codons in your reading frame it will go to the first one. The justification for this is the stop, listed as "*", can't be recognized by the blosum substitution matrix
    inseq2aa=inseq2nt.translate(to_stop=True)

    matrix=matlist.blosum62 #note: depending on the sequence homology it may make sense to use another blosum matrix (or different gap_open gap_close)
    gap_open=-12 #cost to open a gap
    gap_extend=-3 #cost to extend a gap

    alns=pairwise2.align.globalds(inseq1aa, inseq2aa, matrix, gap_open, gap_extend,penalize_end_gaps=0)

    top_aln=alns[0]
    aln_inseq1aa, aln_inseq2aa, score, begin, end = top_aln

#    print 'Aligned Protein Sequences:\n\n'+'Protein Sequence 1 (Reference)\n'+aln_inseq1aa+'\n\n'+'Protein Sequence 2 (to optimize)\n'+aln_inseq2aa+'\n'


    #Module 1.1:Separating inputseqs into codons and adding gaps from the alnAAseq data


    gapsinseq1aa=findOccurences(aln_inseq1aa,"-") #These count the #s of "-" in the global aligned AA sequences
    gapsinseq2aa=findOccurences(aln_inseq2aa,"-")

    ##This loop counts the number of mismatches between amino acid sequences, excluding gaps
    comp_1a_to_2a =zip(aln_inseq1aa,aln_inseq2aa)
    mismatchaacount = 0
    for i,j in comp_1a_to_2a:
        if i!=j and i!='-' and j!='-':
           mismatchaacount +=1
 #   print 'Number of Amino Acid Mismatches: {0}\n'.format(mismatchaacount)


    ##print gapsinseq1aa #this is the locations of the positions of the gaps in each aligned AA sequence
    ##print gapsinseq2aa

    dnagapseq1 = [0]*len(gapsinseq1aa)
    count1= 0
    for i in range(len(dnagapseq1)):
        dnagapseq1[i] = (gapsinseq1aa[i]*3)+count1
        count1-=2

#    print dnagapseq1 #this is the locations of the gaps in the nt sequence

    dnagapseq2 = [0]*len(gapsinseq2aa)
    count2= 0
    for i in range(len(dnagapseq2)):
        dnagapseq2[i]= (gapsinseq2aa[i]*3)+count2
        count2-=2

    ##print dnagapseq2 #see above comment


    rawseq1ntlist=list(rawinseq1nt) #this turns the raw input string nt into a list.
    for x in dnagapseq1:
        rawseq1ntlist.insert(x,'---')
        ''.join(rawseq1ntlist)
    ntgapseq1=''.join(rawseq1ntlist)
    ##print ntgapseq1 #this is the nt string with gaps for each nt of the gapped amino acid (so, one AA gap is expressed as "---")


    rawseq2ntlist=list(rawinseq2nt)
    for x in dnagapseq2:
        rawseq2ntlist.insert(x, '---')
        ''.join(rawseq2ntlist)
    ntgapseq2=''.join(rawseq2ntlist)
    ##print ntgapseq2

    newseq2 = [0]*len(ntgapseq2) #initialize seq to be optimized
    ntgapseq2withchanges = [0]*len(ntgapseq2)

    #Module 2: This section copies codons from seq1 to the new seq for identical amino acids, copies from seq2 for all gaps, and
    #records the locations of all of the mismatched amino acids, and calculates the differences

    #NOTE: the code is currently using the amino acid distances *3 as the mismatch distance, which isn't completely accurate;
    #there is a loop I wrote below that does it more accurately, but I can't get it to ignore mismatches within the same codon, so
    #we currently can't use the more accurate distances... THIS AFFECTS THE OPTIMIZATION ONLY MINIMALLY

    y = 0
    mismatch_AA = [-1] #adds 5' to first gap
    already_ident = 0
    silentmutcount = 0
    ##misandgaps = []

    ##this while loop copies codons from seq1 to newseq for all matches, from seq2 for all gaps, and keeps track of locations of mismatches
    while y < len(aln_inseq1aa):
        if aln_inseq1aa[y] == aln_inseq2aa[y]:  #this copies identicals
            newseq2[(y*3):((y*3)+3)] = ntgapseq1[(y*3):((y*3)+3)]
            ntgapseq2withchanges[(y*3):((y*3)+3)] = ntgapseq1[(y*3):((y*3)+3)]
            if ''.join(map(str, newseq2[(y*3):((y*3)+3)])) != ntgapseq2[(y*3):((y*3)+3)]: #this loop counts what codons actually change
                silentmutcount += 1
            else: already_ident +=1
        elif aln_inseq2aa[y] == '-': #if there's a gap in 2, copy nt of 2
            newseq2[(y*3):((y*3)+3)] = ntgapseq2[(y*3):((y*3)+3)]
            ntgapseq2withchanges[(y*3):((y*3)+3)] = ntgapseq2[(y*3):((y*3)+3)]
    ##        misandgaps.append(y)

        elif aln_inseq1aa[y] == '-': #if there's a gap in 1, copy nt of 2 (redundant, but useful for tracking)
            newseq2[(y*3):((y*3)+3)] = ntgapseq2[(y*3):((y*3)+3)]
            ntgapseq2withchanges[(y*3):((y*3)+3)] = ntgapseq2[(y*3):((y*3)+3)]
    ##        misandgaps.append(y)

        else:
            mismatch_AA.append(y) #makes an int list of the positions of the mismatches in the AA sequences
            ntgapseq2withchanges[(y*3):((y*3)+3)] = ntgapseq2[(y*3):((y*3)+3)]
    ##        misandgaps.append(y)

        y += 1

    newseq2[(y*3):((y*3)+3)] = ntgapseq2[(y*3):((y*3)+3)] #this adds the stop codon at the end from sequence 2

    mismatch_AA.append(len(aln_inseq2aa)-1) #Adds in 3' end for distance calc

    newseqstring = ''.join(map(str, newseq2)) #right now, we print ____ for gaps in seq2 and 000 for mismatches (gaps in seq 1 are just copied codons from seq2)
    ##print newseqstring
    ##print 'mismatched AA positions (with ends):'
    ##print mismatch_AA

    mis_loc = mismatch_AA[1:(len(mismatch_AA)-1)]
    ##print len(mis_loc)

    #misdiffnogap calculates the distance between mismatches in AA seq2, subtracting out the gaps
    misdiffnogap = [((mismatch_AA[n]-mismatch_AA[n-1]-aln_inseq2aa.count('-', mismatch_AA[n-1], mismatch_AA[n]))-1)*3 for n in range(1,len(mismatch_AA))]


    ##print 'mismatches next to each other: {0}'.format(misdiffnogap.count(0))
    ##print 'distance between mismatches, sans gaps:'
    ##print misdiffnogap
    ##print 'number of gaps: {0}'.format(len(misdiffnogap))

    ##loop that adds bases from mismatches that are homologous on either end to the homology lengths
    for n in range(0,len(mis_loc)):
        loc = mis_loc[n]
        triseq1 = ntgapseq1[(loc*3):((loc*3)+3)]
        triseq2 = ntgapseq2[(loc*3):((loc*3)+3)]
        if triseq1[0] == triseq2[0] and triseq1[1] == triseq2[1]:
            misdiffnogap[n] += 2
    ##        print '2L'
        elif triseq1[2] == triseq2[2] and triseq1[1] == triseq2[1]:
            misdiffnogap[n+1] += 2
    ##        print '2R'
        elif triseq1[0] == triseq2[0]:
            misdiffnogap[n] += 1
    ##        print '1L'
        elif triseq1[2] == triseq2[2]:
            misdiffnogap[n+1] += 1
    ##        print '1R'

    num_zeros = misdiffnogap.count(0) #number of amino acids next to each other
    ##print len(misdiffnogap)
    ##print len(mis_loc)


    #Module 3: Copies codons over for single option amino acids (M,W), optimizes amino acids that can change
    #on the right side only, then goes through a hell of an ordeal to optimize amino acids whose codons can change on
    #both sides

    #NOTE: at this time, we did not write any code to try and change the middle of the codon for S, the only
    #codon that could change in the middle sometimes. This occurance will be extremely rare.

    ##Amino Acid Dictionaries
    ##NOTE: Single codon Amino Acids are M,W
    Acid_Dic = dict(
    F = {1:'TTT',2:'TTC'},
    Y = {1:'TAT',2:'TAC'},
    H = {1:'CAT',2:'CAC'},
    Q = {1:'CAA',2:'CAG'},
    N = {1:'AAT',2:'AAC'},
    K = {1:'AAA',2:'AAG'},
    D = {1:'GAT',2:'GAC'},
    E = {1:'GAA',2:'GAG'},
    C = {1:'TGT',2:'TGC'},
    I = {1:'ATT',2:'ATC',3:'ATA'},
    V = {1:'GTT',2:'GTC',3:'GTA',4:'GTG'},
    P = {1:'CCT',2:'CCC',3:'CCA',4:'CCG'},
    T = {1:'ACT',2:'ACC',3:'ACA',4:'ACG'},
    A = {1:'GCT',2:'GCC',3:'GCA',4:'GCG'},
    G = {1:'GGT',2:'GGC',3:'GGA',4:'GGG'},
    R = {1:'CGT',2:'CGC',3:'CGA',4:'CGG',5:'AGA',6:'AGG'},
    S = {1:'TCT',2:'TCC',3:'TCA',4:'TCG',5:'AGT',6:'AGC'},
    L = {1:'TTA',2:'TTG',3:'CTT',4:'CTC',5:'CTA',6:'CTG'})
    R_change = ['F','Y','H','Q','N','K','D','E','C','I','V','P','T','A','G']
    B_change = ['R','S','L']

    #loop for optimization of amino acids
    AA_opt_count_R = 0
    AA_opt_count_L = 0

    tot_aa_change = 0
    tot_aa_copy = 0
    tot_already_opt = 0

    for n in range(0,len(mis_loc)): #for every mismatch location:
        loc = mis_loc[n] #use this instead of n (*3 for DNA location) to call the locations of mismatches
        l_dist = misdiffnogap[n] #left and right distances, respectively
        r_dist = misdiffnogap[n+1]
        AA_to_opt = aln_inseq2aa[loc] #The amino acid in aa sequence2 we need to optimize

        if AA_to_opt == 'M' or AA_to_opt == 'W': #replaces M's and W's, they only have one codon (checked before and after to verify)
            newseq2[(loc*3):((loc*3)+3)] = ntgapseq2[(loc*3):((loc*3)+3)]
            tot_aa_copy += 1
            tot_already_opt += 1

        elif AA_to_opt in R_change: #finds amino acids that can only change on the right
            for s in Acid_Dic[AA_to_opt]:
                poss = Acid_Dic[AA_to_opt][s] #goes through each possible codon in the dictionary for a given AA to optimize
    ##            print "{0},{1}".format(AA_to_opt,s) check to see what aa it's optimizing and how many tries it makes
                if poss[2] == ntgapseq1[((loc*3)+2)]: #looks to see if right end matches the appropriate reference sequence nucleotide
                    newseq2[(loc*3):((loc*3)+3)] = poss #if it does, it switches it
                    if ntgapseq2[(loc*3):((loc*3)+3)] != ''.join(map(str, newseq2[(loc*3):((loc*3)+3)])): #checks to see if the changed amino acid is different than original
                        AA_opt_count_R += 1 #if it does, it adds it to the count of amino acids optimized on the right
                        misdiffnogap[n+1] += 1
                        tot_aa_change += 1
    ##                    print ntgapseq2[(loc*3):((loc*3)+3)] #this would print the old codon and new codon
    ##                    print ''.join(map(str, newseq2[(loc*3):((loc*3)+3)]))
                    else:
                        tot_aa_copy += 1
                        tot_already_opt += 1


        elif AA_to_opt in B_change: #finally, the amino acids that can change on both sides
            for a in Acid_Dic[AA_to_opt]:
                if ''.join(map(str, newseq2[(loc*3):((loc*3)+3)])) == '000': #This stops trying if the codon has already be changed
                    poss = Acid_Dic[AA_to_opt][a]
    ##                print "{0},{1}".format(AA_to_opt,a) check to see what aa it's optimizing and how many tries it makes
                    if poss[2] == ntgapseq1[((loc*3)+2)] and poss[0] == ntgapseq1[((loc*3))]: #can both ends have homology?
                        newseq2[(loc*3):((loc*3)+3)] = poss

                        if ntgapseq2[(loc*3):((loc*3)+3)] != ''.join(map(str, newseq2[(loc*3):((loc*3)+3)])): #is the codon new?
                            tot_aa_change += 1
                            if ntgapseq2[(loc*3)] != ''.join(map(str, newseq2[(loc*3)])): #does it change left side?
                                AA_opt_count_L += 1
                                if ntgapseq2[(loc*3)+1] != ''.join(map(str, newseq2[(loc*3)+1])): #does it change left and middle?
                                    misdiffnogap[n] += 2
                                else:
                                    misdiffnogap[n] += 1

                            if ntgapseq2[(loc*3)+2] != ''.join(map(str, newseq2[(loc*3)+2])): #does it change right side?
                                AA_opt_count_R += 1
                                if ntgapseq2[(loc*3)+1] != ''.join(map(str, newseq2[(loc*3)+1])): #does it change right and middle?
                                    misdiffnogap[n+1] += 2
                                else:
                                    misdiffnogap[n+1] += 1
                        else:
                            tot_aa_copy += 1
                            tot_already_opt += 1

                    else:#these go through in order of most optimal side to change it, then finally at the end says fuck it if neither is optimal and trys to change left then right
                        if abs(l_dist-6) < abs(r_dist-6) and l_dist > 5 and l_dist < 10: # if l_dist is optimal, try to optimize this first
                            if poss[0] == ntgapseq1[((loc*3))]: #left end be fixed?
                                newseq2[(loc*3):((loc*3)+3)] = poss
                                if poss[0] != ntgapseq2[((loc*3))]:
                                    AA_opt_count_L += 1
                                    tot_aa_change += 1
                                else:
                                   tot_already_opt += 1
                                   tot_aa_copy += 1
                        elif abs(l_dist-6) >= abs(r_dist-6) and r_dist > 5 and r_dist < 10: # if r_dist is optimal, try it
                            if poss[2] == ntgapseq1[((loc*3)+2)]: #right end be changed and fixed?
                                newseq2[(loc*3):((loc*3)+3)] = poss
                                if poss[2] != ntgapseq2[((loc*3)+2)]:
                                    AA_opt_count_R += 1
                                    tot_aa_change += 1
                                else:
                                   tot_already_opt += 1
                                   tot_aa_copy += 1
                        elif abs(l_dist-6) < abs(r_dist-6) and l_dist <= 6:
                            if poss[0] == ntgapseq1[((loc*3))]: #left end be changed and fixed?
                                newseq2[(loc*3):((loc*3)+3)] = poss
                                if poss[0] != ntgapseq2[((loc*3))]:
                                    AA_opt_count_L += 1
                                    tot_aa_change += 1
                                else:
                                   tot_already_opt += 1
                                   tot_aa_copy += 1
                        elif abs(l_dist-6) >= abs(r_dist-6) and r_dist <= 6:
                            if poss[2] == ntgapseq1[((loc*3)+2)]: #right end be changed and fixed?
                                newseq2[(loc*3):((loc*3)+3)] = poss
                                if poss[2] != ntgapseq2[((loc*3)+2)]:
                                    AA_opt_count_R += 1
                                    tot_aa_change += 1
                                else:
                                   tot_already_opt += 1
                                   tot_aa_copy += 1
                        else:
                            if poss[0] == ntgapseq1[((loc*3))]: #left end be changed and fixed?
                                newseq2[(loc*3):((loc*3)+3)] = poss
                                if poss[0] != ntgapseq2[((loc*3))]:
                                    AA_opt_count_L += 1
                                    tot_aa_change += 1
                                else:
                                   tot_already_opt += 1
                                   tot_aa_copy += 1
                            elif poss[2] == ntgapseq1[((loc*3)+2)]: #right end be changed and fixed?
                                newseq2[(loc*3):((loc*3)+3)] = poss
                                if poss[2] != ntgapseq2[((loc*3)+2)]:
                                    AA_opt_count_R += 1
                                    tot_aa_change += 1
                                else:
                                   tot_already_opt += 1
                                   tot_aa_copy += 1

        if ''.join(map(str, newseq2[(loc*3):((loc*3)+3)])) == '000':  #This changes the amino acids if they have not been optimized
            newseq2[(loc*3):((loc*3)+3)] = ntgapseq2[(loc*3):((loc*3)+3)]
            tot_aa_copy += 1




    ##print misdiffnogap

    newstring2 = ''.join(map(str, newseq2)) #right now, we print ____ for gaps in seq2 and 000 for mismatches (gaps in seq 1 are just copied codons from seq2)
    ##print newstring2
    for i in range(newseq2.count('-')):
        newseq2.remove('-')
    newstring3 = ''.join(map(str, newseq2))



    comp_org_to_opt =zip(ntgapseq2,newstring2) #this loop compares new seq2 to old and counts changes
    nuc_imp = 0.0
    total_homo = 0.0
    tot_gaps = 0.0
    for i,j in comp_org_to_opt:
        if i!=j:
            nuc_imp += 1

    comp_ref_to_opt =zip(ntgapseq1,newstring2) #this compares homology between ref seq and new sequence
    for i,j in comp_ref_to_opt:
        if i==j and i!='-':
            total_homo += 1
        if i=='-' or j=='-':
            tot_gaps += 1


    z = len(newstring2) - tot_gaps
    perc_improvement = 100.0*(nuc_imp/z)

    perc_homology = 100.0*(total_homo/z)
    print ("\n\nResults:\n{0} codons of matched amino acids in the protein alignment were changed to match reference DNA codons".format(silentmutcount))
    print ("\t {0} matched amino acids already had identical codons\n".format(already_ident))
    print ("{0} codons of mismatched amino acids were changed in total, while {1} were not changed.".format(tot_aa_change,tot_aa_copy))
    print ("\tof those not changed, {0} were already optimal".format(tot_already_opt))
    print ("\t{0} mismatched amino acids in sequence 2 were optimized to add homology to the 3' end".format(AA_opt_count_R))
    print ("\t{0} mismatched amino acids in sequence 2 were optimized to add homology to the 5' end".format(AA_opt_count_L))
    print ("{0} total nucleotides were changed to increase homology between sequences for a total improvement of {1}%, or a total of {2}% nucleotide identities between aligned amino acids".format(nuc_imp,perc_improvement,perc_homology))
    print (newstring3)
    return newstring3

def findOccurences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

