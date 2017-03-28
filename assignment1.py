import mysql.connector
from ast import literal_eval as make_tuple
import os
import pysam
import pybedtools as pbt
import numpy as np
import subprocess

__author__ = 'Josef Moser'


class GoI:
    def __init__(self, name, accession, chrom, start, end, strand, exons, i_start, i_end):
        self.name = name
        self.accession = accession
        self.chrom = chrom[3:]
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = exons
        self.exons_start = str(i_start).lstrip("b'").rstrip(",'").split(",")
        self.exons_end = str(i_end).lstrip("b'").rstrip(",'").split(",")


class Assignment1:
    def __init__(self):
        self.gene = "ATM"
        self.genesfile = "genes.txt"
        self.path = os.path.join(os.getcwd(), self.genesfile)
        self.refseqfile = os.path.join(os.getcwd(),
                                       "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam")
        if not os.path.isfile(self.refseqfile):
            subprocess.call(["wget", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/"
                                     "alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120"
                                     "522.bam"])
        if not os.path.isfile(self.refseqfile + ".bai"):
            subprocess.call(["samtools", "index", self.refseqfile])

        if not os.path.isfile(self.path):
            self.genefh = self.fetch_gene_coordinates("hg19", self.genesfile)

        else:
            self.genefh = open(self.genesfile)

        self.refseq = pysam.AlignmentFile(self.refseqfile, "rb")
        self.goi = self.get_goi(self.genefh)
        self.goi_reads = list(self.refseq.fetch(self.goi.chrom, self.goi.start, self.goi.end))
        self.goi_file = self.write_bam()

    def get_goi(self, fh):
        for line in fh.readlines():
            tup = make_tuple(line)
            if tup[0] == self.gene:
                self.goi = GoI(*[i for i in tup])

        return self.goi

    def write_bam(self):
        fn = self.goi.name + ".bam"
        outfile = pysam.AlignmentFile(fn, "wb", template=self.refseq)
        for read in self.goi_reads:
            outfile.write(read)

        return fn

    @staticmethod
    def fetch_gene_coordinates(genome_reference, file_name):
        print("Connecting to UCSC to fetch data")
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu',
                                      user='genomep',
                                      passwd='password',
                                      db=genome_reference)
        cursor = cnx.cursor()
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        cursor.execute(query)

        fh = open(file_name, "w")
        [fh.write(str(row) + "\n") for row in cursor]
        cursor.close()
        cnx.close()
        print("Done fetching data")

        return fh

    def get_sam_header(self):
        print("\nRefseq Sam Header:")
        for k, v in self.refseq.header["HD"].items():
            if k == "SO":
                print("\tSO (Sorting order of Alignments): {}".format(v))
            if k == "VN":
                print("\tVN (Format version): {}".format(v))
            if k == "GO":
                print("\tGO: (Grouping of alignments): {}".format(v))

    def get_properly_paired_reads_of_gene(self):
        rds = [r for r in self.goi_reads if r.is_paired]  # only paired reads
        print("\nProperly parired reads in gene {}: {}".format(self.goi.accession,
                                                               max([i for i, x in enumerate(rds)])))

    def get_gene_reads_with_indels(self):
        indel_reads = []
        for read in self.goi_reads:
            if not read.is_unmapped:  # if it's mapped
                cigar_line = read.cigar
                for (cigar_type, cigar_length) in cigar_line:
                    if (cigar_type == 1) or (cigar_type == 2):  # insertion or deletion
                        indel_reads.append(read)

        print("\nReads with indels mapping on gene {}: {}".format(self.goi.accession,
                                                                  len(indel_reads)))

    def calculate_total_average_coverage(self):
        # get chromosome length from refseq header:
        g = [{self.goi.chrom: (0, i["LN"])} for i in self.refseq.header["SQ"]
             if i["SN"] == self.goi.chrom][0]
        g = pbt.chromsizes_to_file(g, "genome")
        bt = pbt.BedTool(self.refseqfile)
        bg = bt.genome_coverage(g=g, bg=True)
        cov = np.mean([float(i[3]) for i in bg])
        print("\nThe avarage coverage of Chromosome {} is: {}".format(self.goi.chrom, cov))

    def calculate_gene_average_coverage(self):
        g = {self.goi.chrom: (self.goi.start, self.goi.end)}
        g = pbt.chromsizes_to_file(g, "genome")
        bt = pbt.BedTool(self.goi_file)
        bg = bt.genome_coverage(g=g, bg=True)
        cov = np.mean([float(i[3]) for i in bg])
        print("\nThe avarage coverage of gene {} is: {}".format(self.goi.name, cov))

    def get_number_mapped_reads(self):
        print("\nNumber of mapped reads: {}".format(max([i for i, r in enumerate(self.refseq)
                                                         if not r.is_unmapped])))

    def get_gene_symbol(self):
        print("\nGene symbol for gene {}: {}".format(self.goi.accession, self.goi.name))

    def get_region_of_gene(self):
        print("\nGene region of gene {}: chr{}:{}-{}".format(self.goi.accession, self.goi.chrom,
                                                             self.goi.start, self.goi.end))

    def get_number_of_exons(self):
        print("\nNumber of exons in gene {}: {}".format(self.goi.accession, self.goi.exons))

    def print_summary(self):
        self.get_sam_header()
        self.get_number_mapped_reads()
        self.get_gene_symbol()
        self.get_region_of_gene()
        self.get_number_of_exons()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_gene_average_coverage()
        self.calculate_total_average_coverage()

if __name__ == '__main__':
    print("Assignment 1")
    print(__author__)
    assignment1 = Assignment1()
    assignment1.print_summary()
