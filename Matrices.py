from fractions import Fraction

class Matrices:
    entries = []
    rowsize = 0
    colsize = 0

    def __init__(self, m, n):
        self.setSize(m, n)
        self.entries=[Fraction(i,2) for i in range(1, m*n+1)]

    def setRowSize(self, m):
        self.rowsize = m


    def setColSize(self, n):
        self.colsize = n

    def setSize(self, m, n):
        self.setRowSize(m)
        self.setColSize(n)

    def printRow(self, i):
        print ' '.join(map(str,self.entries[i*self.colsize - self.colsize:i*self.colsize])),

    def printCol(self, j):
        for i in range(1, self.rowsize):
            print self.entries[(i-1)*self.colsize + j - 1]
        print str(self.entries[(self.rowsize-1)*self.colsize + j - 1]),

    def printMatrix(self):
        print '[',
        for i in range(1, self.rowsize):
            self.printRow(i)
            print '\n ',
        self.printRow(self.rowsize),
        print ']'


x=Matrices(3,3)
x.printMatrix()
