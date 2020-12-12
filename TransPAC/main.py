#!/bin/env python
import os
import sys
import pysam
import argparse

from dp.aligner import read_seq_to_aligned_pairs

def non_emptyfile(checkfile):
    return os.path.isfile(checkfile) and os.path.getsize(checkfile) > 0

def check_file_exist(file_):
    if not non_emptyfile(file_):
        msg = "The file %s does not exist" % file_
        raise argparse.ArgumentTypeError(msg)
    return file_

def read_regions(bedfile,flank):
    regions = []
    with open(bedfile,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = max(int(line[1]) - flank,1)
            end = int(line[2]) + flank # might break
            regions.append((chrom,start,end))
    return regions

def skip_read(read,start,end):
    skip = False
    if read.is_secondary:
        skip = True            
    if read.is_supplementary:
        skip = True
    if read.is_unmapped:
        skip = True
    if read.query_sequence == None:
        skip = True
    if read.reference_start > start:
        skip = True
    if read.reference_end < end:
        skip = True
    return skip
        
def extract_sequence(bam,regions,reffn):
    ref = {}
    query = {}
    fasta = pysam.FastaFile(reffn)
    samfile = pysam.AlignmentFile(bam)    
    for chrom,start,end in regions:
        ref_seq = fasta.fetch(reference=chrom,start=start,end=end)
        for index,read in enumerate(samfile.fetch(chrom,start,end)):
            if skip_read(read,start,end):
                continue
            qstart = None
            qend = None
            for ap in read.get_aligned_pairs():
                qpos,rpos = ap
                if rpos is not None and rpos < start and qpos is not None:
                    qstart = qpos
                if rpos is not None and rpos < end and qpos is not None:
                    qend = qpos
                elif rpos is not None and rpos >= end and qpos is not None:
                    pass
            assert qstart is not None
            assert qend is not None
            qseq = read.query_sequence[qstart:qend]
            if len(qseq) == 0:
                continue
            name = "%s_%s_%s_%s" % (index,chrom,start,end)
            ref[name] = ref_seq
            query[name] = qseq
    return (query,ref)

def run_repeatmasker(ins_fasta,outdir):
    command = "RepeatMasker -e rmblast -species human %s -dir %s" % (ins_fasta,outdir)
    os.system(command)

def valid_tsds(table):
    tsds = {}
    with open(table,'w') as fh:
        for line in fh:
            sl = line.strip().split()
            seqid = sl[0]
            tsd1 = sl[2]
            tsd2 = sl[4]
            tsd_ref = sl[6]
            if similar(tsd1,tsd_ref) >0.85 and similar(tsd2,tsd_ref) > 0.85 and len(tsd1)>4 and len(tsd1)<45:
                tsds[seqid] = sl
    return tsds

def ins_seqs(ins_fasta):
    seqs = {}
    ins_seqs = SeqIO.parse(open(ins_fasta),"fasta")
    for i in ins_seqs:
        seqs[i.id]=str(i.seq)
    return seqs
        
def get_mei(outdir,ins_fasta,table):
    meifn = "%s/mei.txt" % outdir
    repeatmasker_output = "%s.out" % ins_fasta
    seqs = ins_seqs(ins_fasta)
    tsds = valid_tsds(table)    
    rmes = parseRepeatMasker(repeatmasker_output,seqs)
    with open(meifn,'w') as fh: 
        for rme in rmes:
            if rme.maybeMEI() is True:
                name = rme.query_id
                if (name in seqs) and (name in tsds):
                    fh.write(rme)
    
def run(args):
    regions = read_regions(args.regions,args.flank)
    query, ref = extract_sequence(args.bam,regions,args.ref)
    alg_output = {
        "table": "%s/table.txt" % args.outdir,
        "ins_fasta": "%s/ins.fasta" % args.outdir
    }
    for fmt in alg_output:
        read_seq_to_aligned_pairs(query,ref,args.k,fmt,alg_output[fmt])
    run_repeatmasker(alg_output["ins_fasta"])
    get_mei(args.outdir,alg_output["ins_fasta"],alg_output["table"])
    
def main():
    parser = argparse.ArgumentParser(description='Detect transduction events')
    parser.add_argument('bam', metavar='bam', type=check_file_exist,
                        help='Alignment file')
    parser.add_argument('ref', metavar='ref', type=check_file_exist,
                        help='Reference fasta file')
    parser.add_argument('regions', metavar='regions', type=check_file_exist,
                        help='BED file with targetted regions')
    parser.add_argument('outdir', metavar='outdir',                        
                        help='output directory')
    parser.add_argument('--flank', metavar="flank", default=1000, type=int,
                        help='Flank length added to regions in BED file')
    parser.add_argument('--k', metavar="k", default=10, type=int,
                        help='Target site duplication length')
    args = parser.parse_args()
    return run(args)
