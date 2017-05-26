import sys
import cmath
import re

# linear map val∈[valmin, valmax] |→ out∈[outmin, outmax]
def scale(val,  valmin,  valmax,  outmin,  outmax):
    return (
        (val - valmin) / (valmax - valmin)
        * (outmax - outmin) + outmin
    )

c         = input("c? (default 0 + 0i):           ") or "0"
crange    = input("c range? (default 1.5):        ") or "1.5"
column    = input("column / row (default 1 1):    ") or "1 1"
cellcount = input("column/row count (default 10): ") or "10"

cellcount = max(int(cellcount), 1)

c = re.sub("\s", "", c)
c = re.sub("i", "j", c)
c = complex(c)

crange = max(float(crange), 0)

inx    = column.find(" ")
row    = int(column[(inx + 1):]) - 0.5
column = int(column[:inx]) - 0.5

outc = (
    c.real
    + scale(column, 0.5, cellcount - 0.5, -crange, crange)
    + (c.imag
    + scale(   row, 0.5, cellcount - 0.5, crange, -crange)
)*1j)

print("\nc = " + str(outc.real) + " + " + str(outc.imag) + "i")
