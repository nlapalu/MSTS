#!/usr/bin/env python

import logging
import os
import sqlite3

class SqliteDB(object):

    def __init__(self, dbfile='', logLevel='DEBUG'):
        
        self.dbfile = dbfile

        if logLevel:
            self.setLoggingLevel(logLevel)

        if dbfile:
            self.conn = sqlite3.connect(dbfile)
        else:
            dbfile='/tmp/tmpSqlite.db' 
            self.conn = sqlite3.connect(dbfile)
            logging.info("SQLiteDB in {}".format(dbfile))

    def setLoggingLevel(self, logLevel):
        logging.basicConfig(level=logLevel)


    def getListOfTables(self):
        pass

    def getDbFileName(self):
        """Return fileName"""

        return self.dbfile

    def  commit(self):
        """Commit transactions"""

        self.conn.commit()


    def deleteDB(self):
        """Delete the DB"""

        self.conn.close()
        os.remove(self.dbfile)
