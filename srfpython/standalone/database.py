import sqlite3, os, glob, sys, traceback
import numpy as np
import pickle
from printcolors import printyellow #, printred


def repartition(z):
    """compute the repartition function of a serie"""
    N = float(len(z))
    zs = np.sort(z)
    x = np.sort(np.unique(zs))
    y = np.ones_like(x)
    i = 0
    for j in xrange(1, len(zs)):
        zz = zs[j]
        if zz == x[i]:
            y[i] += 1.0
        else:
            i += 1
            y[i] = y[i-1] + 1.0
    return x, y / N

########################################################################
def errormsg():
    type, value, trace = sys.exc_info()
    message = "".join(traceback.format_exception(type, value, trace, limit=10))
    return message

########################################################################
class SqliteStdevFunc:
    """user defined aggregate for standard deviation, see database.create_aggregate"""
    def __init__(self):
        self.M  = 0.0
        self.M2 = 0.0
        self.k  = 0
    def step(self, value):
        if value is not None:
            self.M  = (value + self.k * self.M) / (self.k + 1) #mean on the fly
            self.M2 = (value ** 2. + self.k * self.M2) / (self.k + 1) #square mean on the fly
            self.k += 1
    def finalize(self):
        try:
            if self.k > 1:
                if self.M ** 2. > self.M2: return 0. #numerical issues
                return (self.M2 - self.M ** 2.) ** 0.5
            elif self.k == 0:
                return 0.
        except:
            type, value, trace = sys.exc_info()
            message  = "ERROR : user defined aggregate function (%s) failed\n" % self.__class__.__name__
            message += "          "
            message += "          ".join(traceback.format_exception(type, value, trace, limit=5))
            print message
            raise
        return None
########################################################################
class SqliteMedFunc:
    """user defined aggregate for median, see database.create_aggregate"""
    def __init__(self):
        self.X = []
        self.p = 0.5

    def step(self, value):
        if value is not None:
            self.X.append(value)

    def finalize(self):
        try:
            if self.p == 0.5:
                return np.median(self.X)
            if len(self.X) > 3:
                x, y = repartition(self.X)
                return np.interp(self.p, xp = y, fp = x)
            else:
                raise Exception('serie was too short')
            return None
        except:
            type, value, trace = sys.exc_info()
            message  = "ERROR : user defined aggregate function (%s) failed\n" % self.__class__.__name__
            message += "          "
            message += "          ".join(traceback.format_exception(type, value, trace, limit=5))
            print message
            raise
class SqliteP01Func(SqliteMedFunc):
    def __init__(self):
        SqliteMedFunc.__init__(self)
        self.p = 0.01
class SqliteP05Func(SqliteMedFunc):
    def __init__(self):
        SqliteMedFunc.__init__(self)
        self.p = 0.05
class SqliteP16Func(SqliteMedFunc):
    def __init__(self):
        SqliteMedFunc.__init__(self)
        self.p = 0.16
class SqliteP84Func(SqliteMedFunc):
    def __init__(self):
        SqliteMedFunc.__init__(self)
        self.p = 0.84
class SqliteP95Func(SqliteMedFunc):
    def __init__(self):
        SqliteMedFunc.__init__(self)
        self.p = 0.95
class SqliteP99Func(SqliteMedFunc):
    def __init__(self):
        SqliteMedFunc.__init__(self)
        self.p = 0.99

def substring(X, i, j):
    if i == j: return str(X)[i]
    return str(X)[i:j]
########################################################################
########################################################################
#Database will replace database in the future to avoid confusion with module name
class Database(object):
    """
    object to manage sqlite data base,
    data base is connected in the __enter__ method and closed in the __exit__ one,
    use with operator
    """
    #----------------------------------------------------------------
    def __init__(self, sqlitefile, create=False, verbose=True, timeout=120.0):
        """database.__init__ : connect to a sqlite file
            sqlitefile : string, file name to connect to, if not existent (creation), use create = True
            create : bool, set to True if the file needs to be created
            verbose : bool, print messages about the connection
            timeout : float, time after which to give up connexion attempts
        """
        if len(glob.glob(sqlitefile)) != 1 and not create:
            raise Exception('file %s not found, use option create=True to generate a new file' % sqlitefile)
        self.verbose = verbose
        self.sqlitefile = os.path.abspath(sqlitefile)
        self.cnx = None#sqlite3.connect(self.sqlitefile)
        self.transaction = False
        self.lowtransaction = False #for save points
        self.attached = []
        self.timeout = timeout
    #----------------------------------------------------------------
    #MANAGE THE META TABLE (METADATA), objects are stored as pickled strings, use pickle.loads to recover the original objects (see getmeta)
    def dropmeta(self):
        """drop and create the META table"""
        self.cursor.execute('drop table if exists META')
        self.cursor.execute('''create table META (
                KEY      char(100) unique not null,
                VAL      char);''')
    #----------------
    def setmeta(self, key, value):
        """insert an item into the meta table"""
        assert self.transaction
        self.cursor.execute('''
            insert into META (KEY, VAL) values (?, ?)
            ''', (key, pickle.dumps(value)))
    #----------------
    def updatemeta(self, key, value):
        """update an item into the meta table"""
        assert self.transaction
        self.cursor.execute('''
            update META set VAL = ? where KEY = ?
            ''', (pickle.dumps(value), key))
    #----------------
    def getmeta(self, key):
        """load entry from the meta table"""
        cmd = 'select KEY, VAL from META where KEY = ?'
        KEY, VAL = self.cursor.execute(cmd, (key, )).fetchone() 
        return pickle.loads(VAL)
        
    #----------------------------------------------------------------
    def attach(self, sqlitefile, alias):
        """attach another sqlite file to self, see self.detach"""
        self.cursor.execute('''
            attach database "%s" as %s
            ''' % (sqlitefile, alias))
        if self.verbose:
            printyellow("attaching database   : %s" % sqlitefile)
        self.attached.append((sqlitefile, alias))
    #----------------------------------------------------------------
    def disp(self):
        """display the database structure"""
        for w in self.select('select * from sqlite_master'):
            if w[4] is None: continue
            print "#---------------------------"
            print w[4]
    #----------------------------------------------------------------
    def detach(self, alias):
        """detach all attached databases"""
        for bddfile_, alias_ in self.attached:
            if alias == alias_: break
        else: return
        self.cursor.execute('detach database %s' % alias)
        self.attached.remove((bddfile_, alias_))
        if self.verbose:
            printyellow("detaching database   : %s" % bddfile_)
    #----------------------------------------------------------------
    def begintransaction(self):
        """start a new transaction, insersions will be faster, see rollback"""
        if self.transaction:
            raise Exception('a transaction is already open')
        if self.verbose:
            printyellow("starting transaction : %s" % self.sqlitefile)
        self.cursor.execute('begin transaction')
        self.transaction = True
    #----------------------------------------------------------------
    def savepoint(self):
        """leave a save point during the current transaction, rollback may
           remove entries from "begin transaction" or from "last savepoint"
        """
        if not self.transaction:
            raise Exception('no transaction is open')
        if self.verbose:
            printyellow("savepoint            : %s" % self.sqlitefile)
        if self.lowtransaction:
            #release old save point and start a new one
            self.cursor.execute('RELEASE SAVEPOINT LASTSP')            
            self.cursor.execute('SAVEPOINT LASTSP')
        else:
            #start a new savepoint
            self.cursor.execute('SAVEPOINT LASTSP')            
        self.lowtransaction = True
    #----------------------------------------------------------------
    def restarttransaction(self):
        """like savepoint except that uncommited modifications will be commited,
        this might be used if other connections are waiting for their turn to access the database (use timeout >> 1)
        other connection may see the commited changes"""
        if not self.transaction: raise Exception('no transaction is open')
        verbose = self.verbose
        self.verbose = False
        self.commit()
        self.begintransaction()
        self.verbose = verbose
        if self.verbose : 
            printyellow("restart transaction  : %s" % self.sqlitefile)
    #----------------------------------------------------------------
    def rollback(self, crash=True, ignoresavepoint=False):
        """ delete recent insersions
        :param crash: bool, raise an exception after deletions
        :param ignoresavepoint: bool, if True, the rollback includes all entries since begin transaction
        """
        if not self.transaction: 
            raise Exception('no transaction is open')
            
        if self.lowtransaction and not ignoresavepoint:
            if self.verbose:
                printyellow("rolling back         : to last savepoint, %s" % self.sqlitefile)
            self.cnx.execute('''rollback transaction to savepoint LASTSP''')
            self.lowtransaction = False
            self.cnx.commit()
            self.transaction = False            
        else:
            if self.verbose:
                printyellow("rolling back         : to begin transaction, %s" % self.sqlitefile)
            self.cnx.rollback()
            self.transaction = False
        
        if crash:
            printyellow(errormsg())
            raise Exception('transaction failed')

    #----------------------------------------------------------------
    def commit(self):
        """commit entries since begintransaction, closes the current transaction"""
        if self.verbose:
            printyellow("commiting            : %s" % self.sqlitefile)
        if not self.transaction: raise Exception('no transaction in process')
        if self.lowtransaction:
            self.cursor.execute('RELEASE SAVEPOINT LASTSP')            
            self.lowtransaction = False
        self.cnx.commit()
        self.transaction = False
    #----------------------------------------------------------------
    def close(self):
        """"close the database"""
        for sqlitefile, alias in self.attached:
            try: self.detach(alias)
            except : 
                printyellow("could not detach", alias)
                pass  

        if self.transaction: 
            msg = '''you are trying to close a transacting database (%s)
                please chose : 
                1   : commit and close
                2   : rollback since last savepoint and close
                3   : rollback since begin transaction and close''' % self.sqlitefile
            msg = "\n".join([s.strip() for s in msg.split('\n')])
            #choice = generic_nobspy.readinput(msg, timeout = 60., defaultanswer = "2", expectedanswers = ["1", "2", "3"])
            choice = raw_input('%s\n' % msg)
            while not choice in "123":
                choice = raw_input('?')
            if choice == "1":
                self.commit()
            elif choice == "2":
                self.rollback(crash=False, ignoresavepoint=False)
            elif choice == "3":
                self.rollback(crash=False, ignoresavepoint=True)
                
        if self.verbose:
            printyellow("closing              : %s" % self.sqlitefile)
        self.cursor.close()
        self.cnx.close() 
    #---------------------------------------------------
    def __enter__(self):
        """called by the with operator, connect to the database,
        returns self for the "to" keyword"""
        if self.verbose:
            printyellow("connecting to        :", os.path.abspath(self.sqlitefile))
        self.cnx = sqlite3.connect(self.sqlitefile, timeout = self.timeout)
        self.cnx.isolation_level = None
        self.cursor = self.cnx.cursor()
        self.cursor.execute("PRAGMA foreign_keys = ON;") #!!!
        return self
    #---------------------------------------------------    
    def __exit__(self, type, value, trace):
        """whatever happens in the with statement, this will be executed"""
        self.close()     
    #----------------------------------------------------------------
    def insertinto(self, table, dictionaries, insertspec = None):
        """database.insertinto : insert one or more lines into a table
        to insert NULL enter the string "NULL" or None
        """
        if insertspec is not None:
            assert insertspec.lower() in ["or replace", "or rollback", "or abort", "or fail", "or ignore"]
        else:
            insertspec = ""
        if not self.transaction : raise Exception('please run self.begintransaction first')
        if   isinstance(dictionaries, dict): dics = [dictionaries]
        elif isinstance(dictionaries, tuple) or isinstance(dictionaries, list) : dics = dictionaries
        else : raise TypeError('')
        
        for dic in dics:
            keys = []
            vals = []
            for k in dic.keys():
                if dic[k] != 'NULL' and dic[k] != None:
                    keys.append(k)
                    vals.append(dic[k])
            if len(keys) == 0: raise Exception('')
            cmd = "insert %s into %s (" % (insertspec, table)
            for k in keys:
                cmd += '%s, ' % k
            cmd = cmd[:-2] + ')'
            cmd += "values ("    
            cmd += '?, ' * len(keys)
            cmd = cmd[:-2] + ')'
            self.cursor.execute(cmd, tuple(vals))    
            return cmd, vals
    #----------------------------------------------------------------
    def insertorignore(self, *args, **kwargs):
        self.insertinto(insertspec = "or ignore", *args, **kwargs)
    #----------------------------------------------------------------
    def explain(self, cmd, tup=None):  
        "run a selection command in debug mode and yields the query explainations about it"
        cmd = "explain query plan %s" % cmd
        for entry in self.select(cmd, tup = tup):
            yield entry
    #----------------------------------------------------------------
    def select(self, cmd, tup=None, explainquery = False):       
        """database.select : transmit the select command to the database
        return None if the selection is empty
        return a generator otherwise
        example :
            s = mydatabase.select('''select CUSTOMERID, NAME from CUSTOMERS where (NAME = ?)''', ("DUPOND", ))
            if s is None: 
                print "empty selection"
            else: 
                for CUSTOMERID, NAME in s:
                    print CUSTOMERID, NAME
        """
        if explainquery:
            if cmd.split()[0].strip().lower().startswith('select'):
                if tup is not None:
                    for entry in self.cursor.execute('explain query plan %s' % cmd, tup).fetchall():
                        printyellow("eplain query         : %s" % entry[3])
                else:
                    for entry in self.cursor.execute('explain query plan %s' % cmd).fetchall():
                        printyellow("eplain query         : %s" % entry[3])
        #-------
        def generator(selection):
            for item in selection:
                yield item
            selection.close()
            raise StopIteration
        #-------
        
        assert isinstance(cmd, str) or isinstance(cmd, unicode)
        assert tup is None or isinstance(tup, tuple) 
        cmd = cmd.strip()
        if cmd.strip().split()[0].lower() not in ['select', "explain"]:
            raise Exception('cmd must start with key word "select" or "explain", see http://www.sqlite.org/lang_select.html')        
        
        cursortmp = self.cnx.cursor()
        try:
            if tup is not None:
                item0 = cursortmp.execute(cmd, tup).fetchone()
            else:
                item0 = cursortmp.execute(cmd).fetchone()
        except:
            cursortmp.close()
            raise
            
        try:
            if tup is not None:
                selection = cursortmp.execute(cmd, tup)
            else:
                selection = cursortmp.execute(cmd)
        except:
            cursortmp.close()
            raise

        if item0 is None:
            cursortmp.close()
            return None
        else:
            #the generator will close the temporary selection
            #return generator(firstitem=item0, remainingselection=selection)
            return generator(selection)   

    #---------------------------------------------------
    def create_aggregate(self):
        "attach some useful aggregate functions"
        self.cnx.create_aggregate("STD", 1, SqliteStdevFunc) 
        self.cnx.create_aggregate("MED", 1, SqliteMedFunc) 
        self.cnx.create_aggregate("P01", 1, SqliteP01Func) 
        self.cnx.create_aggregate("P05", 1, SqliteP05Func) 
        self.cnx.create_aggregate("P16", 1, SqliteP16Func) 
        self.cnx.create_aggregate("P84", 1, SqliteP84Func) 
        self.cnx.create_aggregate("P95", 1, SqliteP95Func) 
        self.cnx.create_aggregate("P99", 1, SqliteP99Func) 
    #---------------------------------------------------
    def create_function(self):
        "attach some useful functions"
        self.cnx.create_function("LOG", 1, np.log) 
        self.cnx.create_function("SUBSTRING", 3, substring) 
    
if __name__ == "__main__":

    os.system('trash toto.bdd')
    #-----------------
    with Database('toto.bdd', create = True) as bdd:
        bdd.dropmeta()
        bdd.begintransaction()
        try:
            toto = range(10)
            bdd.setmeta('TOTO', toto)
            bdd.commit()
        except:
            bdd.rollback(True)
    #-----------------
    with Database('toto.bdd') as bdd:
        toto = bdd.getmeta('TOTO')
        print toto, type(toto)






