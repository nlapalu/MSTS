#!/usr/bin/env python

class Transcript(object):

    def __init__(self, id, seqid, start, end, strand, gene_id, lCDS=[], lExons=[]):
        """Transcript constructor"""

        self.id = id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id
        self.lCDS = lCDS
        self.lExons = lExons

    def isOnReverseStrand(self):
        """return True if strand -"""

        if self.strand == -1:
            return True
        else:
            return False

    def getLength(self):
        """return total length introns+exons"""

        return self.end-self.start+1

    def getExonTotalLength(self):
        """return the sum of Exon lengths"""

        return sum([exon.end-exon.start+1 for exon in self.lExons])

    def getCDSTotalLength(self):
        """return the length of CDS"""

        return sum([cds.end-cds.start+1 for cds in self.lCDS])

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.id,self.seqid,self.start,self.end,self.strand,self.lCDS, self.lExons,self.gene_id) == (other.id, other.seqid, other.start, other.end, other.strand, other.lCDS, other.lExons,other.gene_id))

    def __repr__(self):
        """Transcript representation"""

        return 'Transcript: {}-{}-{}-{}-{}-{}-{}-{}'.format(self.id,self.seqid,self.start,self.end,self.strand,self.lCDS, self.lExons,self.gene_id)

