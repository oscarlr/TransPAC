#!/usr/bin/python -O
import argparse
import string
import numpy as np
from numpy import array
from Bio import SeqIO
import os
import sys
import aligner_opts

A, C, G, T = 0, 1, 2, 3
int_to_char = {0:'A', 1:'C', 2:'G', 3:'T'}
nucls2num = {'A':0, 'C':1, 'G':2, 'T':3}
neginf = -1000000
indel = -5
scoring = array([[2,-3,-3,-3],
                 [-3,2,-3,-3],
                 [-3,-3,2,-3],
                 [-3,-3,-3,2]])


def vec_pmax_score(a, s):
        "return a + s"
        return a+s

vfunc_pmax_score = np.vectorize(vec_pmax_score)


class AlignmentFinder(object):
        def __init__(self, seq1, seq2):
                self.seq1 = seq1
                self.seq2 = seq2
                self.D = None

        def find_global_alignment(self):
                self.D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                self._compute_array()
                #print(self.D
                return self._traceback()

        def _compute_array(self):
                for i in range(self.seq1.size+1):
                        self.D[i,0] = i*indel
                for j in range(self.seq2.size+1):
                        self.D[0,j] = j*indel
                for i in range(1, self.seq1.size+1):
                        for j in range(1, self.seq2.size+1):
                                self.D[i,j] = max(  self.D[i-1, j-1] + self._get_score(i, j),
                                                    self.D[i-1, j] + indel,self.D[i, j-1] + indel)

        def compute_array(self):
                D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                for i in range(self.seq1.size+1):
                        D[i,0] = i*indel
                for j in range(self.seq2.size+1):
                        D[0,j] = j*indel
                        aligner_opts.compute_array_2dloop(D, scoring, self.seq1.size, self.seq2.size, indel, self.seq1, self.seq2)
                return D

        def compute_array_with_pointers(self):
                '''
              An array Dpointer keeps track of the row you came from
              '''
                D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                Dpointer = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                # alternative would be Dpointer = {}
                for i in range(self.seq1.size+1):
                        D[i,0] = i*indel
                for j in range(self.seq2.size+1):
                        D[0,j] = j*indel
                for i in range(1, self.seq1.size+1):
                        for j in range(1, self.seq2.size+1):
                                D[i,j] = max(  D[i-1, j-1] + self._get_score(i, j),
                                               D[i-1, j] + indel,D[i, j-1] + indel)
                                Dpointer[i,j] = np.argmax(  [D[i-1, j-1] + self._get_score(i, j),
                                                             D[i-1, j] + indel,D[i, j-1] + indel])
                                # alternative would be
                                # if [D[i-1, j-1] + self._get_score(i, j)] == D[i,j]:
                                #      Dpointer[i,j] = (i-1, j-1)
                                # elif D[i-1, j] + indel == D[i,j]:
                                #      Dpointer[i,j] = (i-1, j)
                                # else:
                                #      Dpointer[i,j] = (i, j-1)
                return D, Dpointer



        def _compute_array_mei_helper_tsd(self, i, j, k, Pmax, s):
                #max_v_for_T = neginf # initialize it to a bad value
                #for v in range(j, min(j+k, self.seq2.size+1)): # need to ensure k doesn't exceed the max length of seq2
                #		curr_v = Pmax[i-1, v] + self._get_score(i,j)
                #		if curr_v > max_v_for_T:
                #max_v_for_T = curr_v
                #s = self._get_score(i,j)
                #max_v_for_T = max(map(lambda v: Pmax[i-1, v] + s, range(j, min(j+k, self.seq2.size+1)))) # need to ensure k doesn't exceed the max length of seq2
                return np.max(vec_pmax_score(Pmax[i-1, j:(j+k)], s))


        #vfunc_pmax_score = np.vectorize(vec_pmax_score)

        #	return max_v_for_T

        def _compute_array_mei(self, k):
                '''
              Need to fill out the matrices neccessary to detect and insertion AND
              a TSD duplication (note the dfam hmms are still be ing used to detect
              the mobile elements).
              k = maximal size of our tsd

              Matrix definitions:
              J = jump penalty
              P = Standard dynamic programming (use the _compute_array function above)
              Pmax = the max value at column j up to point i (prefix)
              Pmax(i, j) = max(P(i,j), Pmax(i-1,j)
              T = Looks up k-bases in Dmax (e.g. k columns forward) at the position i-1
              for best jump (TSD)
              T(i, j) = max_{j < v < j+k, u < i} (P(u,v), + scoring(seq1[i-1], seq2[j-1]) - J)
              for speed we write this as:
              T(i, j) = max_{j < v < j+k} (Pmax(i-1, v), + scoring(seq1[i-1], seq2[j-1]) - J)		
              initialization condition
              T(0,j) = 0
              S = allows a jump from T and standard DP options (suffix)
              S = max(T(i-1,j-1) + scoring(seq1[i-1], seq2[j-1]), .... <standard DP operations)
              '''
                P = self.compute_array()
                Pmax =  np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                T =  np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                S =  np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                # pointers to keep track of stufff
                #Tpointer = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                #	Spointer = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int_)
                # fill out Pmax, T
                Pmax[0,0] = P[0,0]
                S[0,0]    = neginf
                T[0,0]    = neginf  # we might want to change this to negative infinity
                for j in range(1, self.seq2.size+1):
                        Pmax[0,j] = P[0,j]
                        T[0,j]    = neginf
                        S[0,j]    = neginf
                for i in range(1, self.seq1.size+1):
                        Pmax[i,0] = max(P[i,0], Pmax[i-1,0])
                        # find maximum value in previous Pmax row within k
                        T[i,0]    = neginf
                        S[i,0]    = neginf
                '''
              for i in range(1, self.seq1.size+1):
                      for j in range(1, self.seq2.size+1):
                              Pmax[i,j] = max(P[i,j], Pmax[i-1,j])
                              score_ij  = self._get_score(i,j)
                              # for T find maximum value in previous Pmax row within k				
                              T[i,j]    = self._compute_array_mei_helper_tsd(i, j, k, Pmax, score_ij)
                              S[i,j] = max(T[i-1,j-1]+score_ij, max(  S[i-1, j-1] + score_ij,
                                                  S[i-1, j] + indel,S[i, j-1] + indel))
              '''
                aligner_opts.compute_array_mei_2dloop (P, Pmax, T, S, scoring, self.seq1.size, self.seq2.size, indel, self.seq1, self.seq2, k)

                return  P, Pmax, T, S

        def _get_score(self, i, j):
                ''' The indexing is quite tricky because the matrix as one more row & column.
              That causes a shift between the matrix index and the sequence indices.
              Therefore, to obtain the correct nucleotide in the sequence, we must
              substract 1 to the matrix index. '''
                return scoring[self.seq1[i-1], self.seq2[j-1]]
        def _get_aligned_pair(self, i, j):
                n1 = int_to_char[self.seq1[i-1]] if i>0 else '_'
                n2 = int_to_char[self.seq2[j-1]] if j>0 else '_'
                return (n1, n2)
        def _traceback(self):
                alignment= []
                i = self.seq1.size
                j = self.seq2.size
                while i >0 and j>0:
                        if self.D[i-1, j-1] + self._get_score(i, j) == self.D[i,j]:
                                alignment.append(self._get_aligned_pair(i, j))
                                i -= 1
                                j -= 1
                        elif self.D[i-1, j] + indel == self.D[i,j]:
                                alignment.append(self._get_aligned_pair(i, 0))
                                i -= 1
                        else:
                                alignment.append(self._get_aligned_pair(0, j))
                                j -= 1
                while i > 0:
                        alignment.append(self._get_aligned_pair(i, 0))
                        i -= 1
                while j > 0:
                        alignment.append(self._get_aligned_pair(0, j))
                        j -= 1
                alignment.reverse()
                return alignment  

        def _traceback_mei(self, P, Pmax, T, S, k):
                alignment = []
                states    = [] # matrices that you are in 
                features  = [] # sequence you are in (Prefix (P), TSD (T), Insertion(I), TSD (T), Suffix (S))
                #start off from the end of the suffix array
                i = self.seq1.size
                j = self.seq2.size
                # start S traceback
                while i >0 and j>0:
                        # check to see if it came from the T matrix
                        if T[i-1,j-1] + self._get_score(i,j) == S[i,j]:
                                alignment.append(self._get_aligned_pair(i, j))
                                states.append("S")
                                features.append("S")
                                i -= 1
                                j -= 1
                                break # exit from the loop
                        elif S[i-1, j-1] + self._get_score(i, j) == S[i,j]:
                                alignment.append(self._get_aligned_pair(i, j))
                                states.append("S")
                                features.append("S")
                                i -= 1
                                j -= 1
                        elif S[i-1, j] + indel == S[i,j]:
                                alignment.append(self._get_aligned_pair(i, 0))
                                states.append("S")
                                features.append("S")
                                i -= 1
                        else :
                                alignment.append(self._get_aligned_pair(0, j))
                                states.append("S")
                                features.append("S")
                                j -= 1

                # now traverse the jump array
                alignment.append(self._get_aligned_pair(i, j))
                states.append("T")
                features.append("T")
                tsd_size = 1 # initialize tsd
                for v in range(j,min(j+k, self.seq2.size+1)):	
                        # identify target site duplication (TSD)
                        if Pmax[i-1, v] + self._get_score(i,j) == T[i,j]:
                                tsd_size = (v-j)+1 # v is where T jumps to S 
                                for t in range(tsd_size):
                                        features[-t-1] = "T" # change annotation of suffix to a TSD
                                j  = v
                                # now find how far back in P we need to jump
                                u  = i-1
                                while (P[u,v] != Pmax[i-1,v]): 
                                        # put in a gap while we are in an MEI
                                        n1 = int_to_char[self.seq1[u]]
                                        n2 = "_"
                                        alignment.append((n1,n2))
                                        states.append("T")
                                        features.append("I")
                                        u -= 1
                                i = u # set i = u to continue for P

                                break # exit for the loop and move to prefix			
                # now traverse the initial P array
                while i >0 and j>0:
                        if P[i-1, j-1] + self._get_score(i, j) == P[i,j]:
                                alignment.append(self._get_aligned_pair(i, j))
                                states.append("P")
                                if tsd_size != 0: # I've filled in the TSD length
                                        features.append("T")
                                        tsd_size -= 1
                                else:
                                        features.append("P")
                                i -= 1
                                j -= 1
                        elif P[i-1, j] + indel == P[i,j]:
                                alignment.append(self._get_aligned_pair(i, 0))
                                states.append("P")
                                features.append("P")
                                i -= 1
                        else:
                                alignment.append(self._get_aligned_pair(0, j))
                                states.append("P")
                                features.append("P")
                                j -= 1
                while i > 0:
                        alignment.append(self._get_aligned_pair(i, 0))
                        states.append("P")
                        features.append("P")
                        i -= 1
                while j > 0:
                        alignment.append(self._get_aligned_pair(0, j))
                        states.append("P")
                        features.append("P")
                        j -= 1
                alignment.reverse()
                states.reverse()
                features.reverse()
                return alignment, states, features

def print_sequences(pairs):
        top_seq = []
        bottom_seq = []
        for (t,b) in pairs:
                bottom_seq.append(b)
                top_seq.append(t)
        for n in top_seq:
                print(n,)
        print(' ')
        for n in bottom_seq:
                print(n,)
        print("")


def get_seq_from_pairs (pairs, ind=0):
        '''
      ind correspond to 0 for query and 1 for ref
      '''
        query_list = [x[ind] for x in pairs]
        return "".join(filter(lambda x: x != "_", query_list))


def read_seq_to_aligned_pairs (fseq, cseq, k, table_fmt, ins_fmt, insertion_size=10000):
        # read_seq_to_aligned_pairs(query,ref,args.k,fmt,alg_output[fmt])
        # full = SeqIO.parse(open(fn1),"fasta")
        # combine = SeqIO.parse(open(fn2),"fasta")
        # fseq = {}
        # cseq ={}
        # for f in full:
        # 	fseq[f.id] = str(f.seq)

        # for c in combine:
        # 	cseq[c.id] = str(c.seq)

        ins_outputfh = open(ins_fmt,'w')
        table_outputfh = open(table_fmt,'w')

        table_outputfh.write("seqid\ttsd1_len\ttsd1\ttsd2_len\ttsd2\ttsdref_len\ttsdref\tins_len\tinsertion\tscore\n")

        for seqid in fseq:
                print("running %s. LEN=%s ..." % (seqid,len(fseq[seqid])),file=sys.stdout)
                fnucs = fseq[seqid].upper()
                cnucs = cseq[seqid].upper()
                if "N" in fnucs or "N" in cnucs:
                        #print(>>sys.stderr, "...'N' found so ignoring %s" %(seqid))
                        continue
                nucls_num = list(map(lambda x: nucls2num[x],fnucs))
                #print(nucls_num)
                insertion = array(nucls_num, dtype=np.int16)
                ref = array(list(map(lambda x: nucls2num[x], cnucs)), dtype=np.int16)
                if len(insertion) > insertion_size:
                        continue
                aligner_test = AlignmentFinder(insertion, ref)
                #pairs_test = aligner_test.find_global_alignment()
                #print(>>sys.stderr, "running mei alignment ...")
                P, Pmax, T, S = aligner_test._compute_array_mei (k)
                #print(>>sys.stderr, "printing mei features ...")
                pairs_test, states_test, features_test = aligner_test._traceback_mei(P, Pmax, T, S, k)

                f_str = "".join(features_test)	
                # find the start and stop index in the alignment 
                min_index = f_str.find("I")
                max_index = f_str.rfind("I")
                # return the subarray [(q1, r1), (q2, r2)...(q_n, r_n)]
                sub_pairs_list = pairs_test[min_index:max_index+1]
                ins_seq = get_seq_from_pairs (sub_pairs_list)
                ins_len = len(ins_seq)

                if len(ins_seq) > 0:
                        ins_outputfh.write(">%s\n%s\n" %(seqid, ins_seq))
                        
                t1s, t1e = -1,-1
                t2s, t2e = -1,-1
                for i in range(len(features_test)):
                        if features_test[i] == "T": # check where T is
                                if t1s == -1: # if you haven't set the first TSD - set
                                        t1s,t1e = i,i+1 
                                elif t2s == -1: # if you haven't set the second
                                        if t1e == i: # check if you're still in the first
                                                t1e = i+1
                                        else: # otherwise set the second
                                                t2s,t2e = i,i+1
                                else: # otherwise you are in the second already
                                        t2e = i+1
                tsd1 = get_seq_from_pairs (pairs_test[t1s:t1e])
                tsd1_len = len(tsd1)
                tsd2 = get_seq_from_pairs (pairs_test[t1s:t1e])
                tsd2_len = len(tsd2)
                tsdref = get_seq_from_pairs (pairs_test[t1s:t1e], 1)
                tsdref_len = len(tsdref)
                table_outputfh.write("%s\n" % "\t".join(map(str, [seqid,tsd1_len, tsd1, tsd2_len, tsd2,tsdref_len, tsdref,ins_len, ins_seq, S[-1][-1]])))

        ins_outputfh.close()
        table_outputfh.close()
                # else:
                # 	print(len(insertion)
                # 	print(len(ref)
                # 	print_sequences(pairs_test)
                # 	print(" ".join(states_test)
                # 	print(" ".join(features_test)	


# if __name__ == "__main__":
#     if len(sys.argv) > 2:
# 	    parser = argparse.ArgumentParser(description='Simple Dynamic Programming Aligner for performing MEI with TSD')

# 	    parser.add_argument('fn1', type=str,
#                    help='query file in fasta format')
# 	    parser.add_argument('fn2', type=str,
# 		    help='ref file in fasta format')
# 	    parser.add_argument('k', type=int,
# 		    help='k (i.e. max size allowed for TSD)')
# 	    parser.add_argument('--fmt', default='full', help='output format')
# 	    args = parser.parse_args()

# 	    fn1 = args.fn1
# 	    fn2 = args.fn2
# 	    k = int(args.k)
# 	    read_file_to_aligned_pairs(fn1, fn2, k, args.fmt)
# 	    sys.exit()

#     s1 = array([G, T, A, C, A, G, T, A], dtype=np.int16)
#     s2 = array([G, G, T, A, C, G, T], dtype=np.int16)
#     print("Initial Sequences"
#     print("".join(map(lambda x: int_to_char[x], s1))
#     print("".join(map(lambda x: int_to_char[x], s2))


#     aligner = AlignmentFinder(s1, s2)
#     pairs = aligner.find_global_alignment()
#     print_sequences(pairs)
#     # try an mei detection
#     k = 5
#     prefix="GC"
#     tsd="TACC"
#     insertion="GATCA"
#     suffix="TTTAT"
#     s1 = map(lambda x: nucls2num[x], prefix+tsd+insertion+tsd+suffix)
#     s2 = map(lambda x: nucls2num[x], prefix+tsd+suffix)
#     s1 = array(s1, dtype=np.int16)
#     s2 = array(s2, dtype=np.int16)
#     print("tsd: ", tsd
#     #s1 = array([G, C, T, A, C, C, G, A, T, C, A, T, A, C, C, T, T, T, A, T], dtype=np.int16)
#     #s2 = array([G, C, T, A, C, C, T, T, T, A, T], dtype=np.int16)

#     print("Initial Sequences"
#     print("".join(map(lambda x: int_to_char[x], s1))
#     print("".join(map(lambda x: int_to_char[x], s2))
#     aligner2 = AlignmentFinder(s1, s2)
#     pairs2a = aligner2.find_global_alignment()
#     P, Pmax, T, S = aligner2._compute_array_mei (k)
#     print(P
#     #print(T
#     #print(S
#     pairs2, states2, features2 = aligner2._traceback_mei(P, Pmax, T, S, k)
#     print_sequences(pairs2)
#     print(" ".join(states2)
#     print(" ".join(features2)
#     print("Score: ",  S[len(aligner2.seq1), len(aligner2.seq2)]

#     #print_sequences(pairs)
#     # now print(just P
#     print("printing alignment of just the P array"
#     print_sequences(pairs2a)


