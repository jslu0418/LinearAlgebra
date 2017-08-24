###########################################################################
# This script is used to implement essential gaussianElimination          #
# Will be improved to a new edition with gaussian-jordan elimination.     #
# (Carl Friedrich Gauss and Wilhelm Jordan)                               #
# Possible bugs: haven't done fraction number process when you define     #
# your own matrix at line, and haven't done enough boundary verification  #
###########################################################################
#############################
# Author: Jiashuai Lu       #
# Email: jslu0418@gmail.com #
#############################

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

    def gaussianElimination(self):
        for i in range(1, self.rowsize + 1): # for all rows in M #
            if i == 1: # if it's first row #
                p = self.nonzeroPositionInRow(1)
                if p == 1: # if nonzeropos equal 1 go to next row #
                    self.gaussianEliminationMultiply(1)
                    continue
                else:
                    if self.rowsize == 1: # if only one row, stop #
                        self.gaussianEliminationMultiply(1)
                        break
                    else:
                    # find the row with nonzeroposion at 1 or  #
                    #   the earliest nonzeropostion row at all #
                        earliestPos = p if p is not None else self.colsize+1
                        earliestPosRow = 1
                        for i2 in range(2, self.rowsize + 1):
                            q = self.nonzeroPositionInRow(i2)
                            if q == 1:
                                self.gaussianEliminationExchange(1, i2)
                                break
                            else:
                                if q > foremost:
                                    earliestPos = q
                                    earliestPosRow = i2
                        if earliestPosRow != 1:
                            self.gaussianEliminationExchange(1, earliestPosRow)
                        self.gaussianEliminationMultiply(1)
            else:
                #######################################################################
                # not the first row, there are 3 possible situation                   #
                # 1. this row has nonzeroposition exactly equal to last row's nzp + 1 #
                # 2. after onestep elimination this row accord to the rule 1          #
                # 3. if not 2, check following rows find a row with nzp equal to last #
                #     row's nzp + 1 or with the least nzp at all following rows.      #
                #######################################################################
                p = self.nonzeroPositionInRow(i-1)
                if p == self.colsize or p == None:
                    break
                earliestPos = self.colsize + 1
                earliestPosRow = None
                for i3 in range(i, self.rowsize + 1):
                    q = self.nonzeroPositionInRow(i3)
                    if q == p + 1:
                        earliestPos = q
                        earliestPosRow = i3
                        break
                    else:
                        self.gaussianEliminationOneStep(i3)
                        q = self.nonzeroPositionInRow(i3)
                        if q == p + 1:
                            earliestPos = q
                            earliestPosRow = i3
                            break
                        else:
                            if q is not None and q < earliestPos:
                                earliestPos = q
                                earliestPosRow = i3
                if earliestPosRow is not None and earliestPosRow != i:
                    self.gaussianEliminationExchange(i3, earliestPosRow)
            self.printMatrix()
        return True

    def nonzeroPositionInRow(self, i):
        count=0
        for e in self.entries[i*self.colsize - self.colsize:i*self.colsize]:
            count = count + 1
            if e != 0:
                return count
        return None

    def gaussianEliminationOneStep(self, i):
        for j in range(1, self.colsize + 1):
            if self.entries[i*self.colsize - self.colsize + j - 1] != 0:
                for k in range(1, i): # k's most value should be h #
                    nzp = self.nonzeroPositionInRow(k)
                    if nzp is not None and nzp == j:
                        if self.entries[k*self.colsize -self.colsize + nzp -1] == 1:
                            self.gaussianEliminationAdd(i, k)
                            self.gaussianEliminationMultiply(i)
        return True

    def gaussianEliminationExchange(self, i, j):
        temp = self.entries[i*self.colsize - self.colsize:i*self.colsize]
        for k in range(1, self.rowsize + 1):
            self.entries[i*self.colsize - self.colsize + k-1] = self.entries[j*self.colsize - self.colsize + k-1]
            self.entries[j*self.colsize - self.colsize + k-1] = temp[k-1]

    def gaussianEliminationAdd(self, i, j):
        nzp = self.nonzeroPositionInRow(j)
        mlp = self.entries[i*self.colsize - self.colsize + nzp - 1] / self.entries[j*self.colsize-self.colsize + nzp -1]
        for k in range(0, self.colsize):
            self.entries[i*self.colsize - self.colsize + k] = self.entries[i*self.colsize - self.colsize + k] - self.entries[j*self.colsize - self.colsize + k] * mlp
        return True

    def gaussianEliminationMultiply(self, i):
        nzp = self.nonzeroPositionInRow(i)
        if nzp is not None:
            divisor = self.entries[i*self.colsize - self.colsize + nzp - 1]
            for k in range(0, self.colsize):
                self.entries[i*self.colsize - self.colsize + k] = self.entries[i*self.colsize - self.colsize + k] / divisor
        return True

x=Matrices(3,4)
x.entries=[1, -2, 3, 9, -1, 3, 0, -4, 2, -5, 5, 17]
x.printMatrix()
x.gaussianElimination()
x.printMatrix()
