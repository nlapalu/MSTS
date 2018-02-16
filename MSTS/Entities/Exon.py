#!/usrbin/env python

class Exon(object):


    def __init__(self,exon_id, seqid, start, end, strand, lTranscript_ids):
        """Exon constructor"""

        self.exon_id = exon_id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.lTranscript_ids = lTranscript_ids

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.exon_id,self.seqid,self.start,self.end,self.strand,self.lTranscript_ids) == (other.exon_id,other.seqid, other.start, other.end, other.strand, other.lTranscript_ids))

    def __repr__(self):
        """Exon representation"""

        return 'Exon: {}-{}-{}-{}-{}-{}'.format(self.exon_id,self.seqid,self.start,self.end,self.strand,self.lTranscript_ids)


