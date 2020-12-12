#cython: profile=True
# cython: wraparound=False
# cython: boundscheck=False
# cython: cdivision=True

import numpy as np
cimport numpy as np

ctypedef np.float_t DTYPE_FLOAT_t
ctypedef np.int_t DTYPE_t
ctypedef np.int16_t DTYPE_int16
ctypedef np.uint16_t dtype_tuint

cdef inline int int_min(int a, int b): return a if a <= b else b

def compute_array_2dloop (np.ndarray[DTYPE_t, ndim=2] D, np.ndarray[DTYPE_t, ndim=2] s, int n1, int n2, int indel, np.ndarray[DTYPE_int16] seq1, np.ndarray[DTYPE_int16] seq2):
    cdef int i = 1
    cdef int j = 1
    cdef int match   = 0
    cdef int ins_gap = 0
    cdef int del_gap = 0
    for i from 1 <= i < n1+1:
        for j from 1 <= j < n2+1:
            match     = D[i-1, j-1] + s[seq1[i-1],seq2[j-1]]
            ins_gap   = D[i-1, j] + indel
            del_gap   = D[i, j-1] + indel
            if match > ins_gap:
                if match > del_gap:
                    D[i,j] = match
                else:
                    D[i,j] = del_gap
            elif ins_gap > del_gap:
                D[i,j] = ins_gap
            else:
                D[i,j] = del_gap


def compute_array_mei_2dloop (np.ndarray[DTYPE_t, ndim=2] P, np.ndarray[DTYPE_t, ndim=2] Pmax, np.ndarray[DTYPE_t, ndim=2] T, np.ndarray[DTYPE_t, ndim=2] S, np.ndarray[DTYPE_t, ndim=2] scoring, int n1, int n2, int indel, np.ndarray[DTYPE_int16] seq1, np.ndarray[DTYPE_int16] seq2, int k):
     cdef int i = 1
     cdef int j = 1
     cdef int match   = 0
     cdef int ins_gap = 0
     cdef int del_gap = 0
     cdef int score_ij = 0
     cdef int curr_v = 0
     cdef int max_v_for_T = 0
     cdef int Smax = 0
     cdef int vmax = 0
     cdef int prev_max_v_for_T = 0
     cdef int prev_max_v_for_T_index = 0 
     for i from 1 <= i < n1+1:
         for j from 1 <= j < n2+1:
             Pmax[i,j] = max(P[i,j], Pmax[i-1,j])
             score_ij  = scoring[seq1[i-1],seq2[j-1]]
             # for T find maximum value in previous Pmax row within k			 
             # first check if you need to update k
             vmax = min(j+k, n2+1)
             if prev_max_v_for_T_index >= j:
                 max_v_for_T = Pmax[i-1, vmax-1]
                 if prev_max_v_for_T > max_v_for_T:
                     max_v_for_T = prev_max_v_for_T
                 else:
                     prev_max_v_for_T = max_v_for_T
                     prev_max_v_for_T_index = vmax-1
             else:
                 max_v_for_T = Pmax[i-1, j]
                 for v from j+1 <= v < vmax:
                     curr_v = Pmax[i-1, v]
                     if curr_v > max_v_for_T:
                         max_v_for_T = curr_v
                         prev_max_v_for_T_index = v
                 prev_max_v_for_T = max_v_for_T
             # set T[i,j]
             T[i,j]    = max_v_for_T + score_ij
             # compute S
             match     = S[i-1, j-1] + score_ij
             ins_gap   = S[i-1, j] + indel
             del_gap   = S[i, j-1] + indel
             if match > ins_gap:
                if match > del_gap:
                    S_max = match
                else:
                    S_max = del_gap
             elif ins_gap > del_gap:
                 S_max = ins_gap
             else:
                 S_max = del_gap

             S[i,j] = max(T[i-1,j-1]+score_ij, S_max)
                

