#!/usrbin/env python

class Feature(object):

    def __init__(self, id, seqid, source, type, start, end, score, strand, phase, attributes):
        """Feature constructor"""

        self.id = id
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __eq__(self, other):
        """Equality on all args"""
      
        return ((self.id,self.seqid,self.source,self.type,self.start,self.end,self.score,self.strand,self.phase,self.attributes) == (other.id, other.seqid, other.source,other.type, other.start, other.end, other.score,other.strand, other.phase,other.attributes))

    def __repr__(self):
        """Feature representation"""

        return 'Feature: {}-{}-{}-{}-{}-{}-{}-{}-{}-{}'.format(self.type,self.id,self.seqid,self.source,self.start,self.end,self.score, self.strand,self.phase,self.attributes)
