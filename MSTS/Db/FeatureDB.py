#!/usr/bin/env python

import os
import logging
import sqlite3
from pickle import dumps, loads

from MSTS.Db.SqliteDB import SqliteDB

from MSTS.Entities.Feature import Feature


class FeatureDB(SqliteDB):

    def __init__(self, dbfile='', noCreate=False,logLevel='ERROR'):
        """FeatureDB Constructor"""

        SqliteDB.__init__(self, dbfile=dbfile, logLevel=logLevel)
        if not noCreate:
            self._createDBSchema()

    def _createDBSchema(self):
        """Create Database Schema"""

        self._createFeatureTable()
        self._createFeatureIndexes()

    def _createFeatureTable(self):
        """Create table Feature"""

        self.conn.execute('''create table feature (
                                                   id text not null,
                                                   seqid text not null,
                                                   source text not null,
                                                   type text not null,
                                                   start int not null,
                                                   end int not null,
                                                   score real,
                                                   strand int, 
                                                   phase int,
                                                   attributes text);''')
        self.conn.commit()

    def _createFeatureIndexes(self):
        """Create Indexes in Feature table"""

        self.conn.execute('''create index seqid_gene_idx on feature(seqid,type,start,end);''')
        self.conn.commit()


    def insertlFeatures(self, lFeatures):
        """Insert a list of features"""

        logging.info('Inserting features in db')

        self.conn.executemany('''insert into feature(id,seqid,source,type,start,end,score,strand,phase,attributes) values (?,?,?,?,?,?,?,?,?,?)''', [(f.id,f.seqid,f.source,f.type,f.start,f.end,f.score,f.strand,f.phase,dumps(f.attributes)) for f in lFeatures])
        self.conn.commit()

        logging.info('{} feature(s) inserted in db'.format(len(lFeatures)))

        return len(lFeatures)


    def selectAllFeatures(self):
        """Select All Features"""

        lFeatures = []
        cursor = self.conn.execute('''select id,seqid,source,type,start,end,score,strand,phase,attributes from feature''')
        for row in cursor:
            lFeatures.append(Feature(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],loads(row[9])))

        return lFeatures


    def selectFeatureTypeFromReference(self, chrom, type):
        """Select features from a specific reference"""

        lFeatures = []
        cursor = self.conn.execute('''select id,seqid,source,type,start,end,score,strand,phase,attributes from feature where seqid = \'{}\' and type = \'{}\' '''.format(chrom, type))
        for row in cursor:
            lFeatures.append(Feature(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],loads(row[9])))

        return lFeatures


    def selectFeatureTypeFromCoordinates(self, type, chrom, start, end):
        """Select features from specific references"""

        lFeatures = []
        cursor = self.conn.execute('''select id,seqid,source,type,start,end,score,strand,phase,attributes from feature where seqid = \'{}\' and type = \'{}\' and ((end > {} and end < {}) or (start < {} and end > {})) order by start'''.format(chrom, type, start, end, start, end))
        for row in cursor:
            lFeatures.append(Feature(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],loads(row[9])))

        return lFeatures

    def selectFeatureFromIdListAndType(self, chrom, lIds, featType):
        """Select features from Id"""

        lFeatures = []
        cursor = self.conn.execute('''select id,seqid,source,type,start,end,score,strand,phase,attributes from feature where seqid = \'{}\' and id in ({}) and type = \'{}\''''.format(chrom, ",".join(['"{}"'.format(x) for x in lIds]),featType))
        for row in cursor:
            lFeatures.append(Feature(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],loads(row[9])))

        return lFeatures

    def selectReferences(self):
        """Select All References"""

        lReferences = []
        cursor = self.conn.execute('''select distinct seqid from feature''')
        for row in cursor:
            lReferences.append(row[0])

        return lReferences


