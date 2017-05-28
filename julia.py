from datetime import datetime
import time
import math
import cmath
import re
import codecs
import sys

# why arent these in math in the first place?
def sec(x):
    return 1 / math.cos(x)

def csc(x):
    return 1 / math.sin(x)

def cot(x):
    return 1 / math.tan(x)

def sign(x):
    return math.copysign(1, x)

math.sec = sec
math.csc = csc
math.cot = cot
math.sign = sign

def eval_fn(z, c):
    global fncode
    ret: float
    try:
        ret = eval(fncode)
    except (ValueError, OverflowError, ZeroDivisionError):
        # negative number in a log or root probably
        ret = float("nan")
    return ret

# linear map val∈[valmin, valmax] |→ out∈[outmin, outmax]
def scale(val,  valmin,  valmax,  outmin,  outmax):
    return (
        (val - valmin) / (valmax - valmin)
        * (outmax - outmin) + outmin
    )

def stgX(x):
    global graph, columnwidth
    return scale(x, 0, columnwidth, graph['x']['min'], graph['x']['max'])

def stgY(y):
    global graph, rowheight
    return scale(y, 0, rowheight, graph['y']['max'], graph['y']['min'])

def clamp(v, lo, hi):
    return max(min(v, hi), lo)

def signstr(num):
    if math.sign(num) == -1:
        return "-"
    else:
        return "+"

print("enter a formula to plot in terms of z and c ({z, c} ∈ ℂ)")

def get_fn():
    return input("z_n+1(z, c) := ")

orig_fn = fn = get_fn() or "z**2 + c"

# replace stuff like 2tan(4x) with 2*tan(4*x)
fn = re.sub(r"(\d+)([a-zA-Z]+)", r"\1 * \2", fn)

# ln = log
fn = re.sub(r"ln", r"log", fn)

# when you type (x + 2)(x - 2) you probably meant to multiply them right?
fn = re.sub(r"\)\s*\(", r")*(", fn)

fn = re.sub(r"π", r"pi", fn)

fn = re.sub(r"(?<!cmath\.)\b(pi|e|tau|inf|infj|nan|nanj)\b", r"cmath.\1", fn)

# (?<! ...)

# sinz, sin z, cos c, etc.
fn = re.sub(r"""(?<!cmath\.)\b(phase|polar|exp|log10|sqrt|acos|asin|atan|cos|
|sin|tan|acosh|asinh|atanh|cosh|sinh|tanh|isfinite|isinf|isnan|log|
|rect)\s*([zc])\b""", r"\1(\2)", fn)

# sin(z) ...
fn = re.sub(r"""(?<!cmath\.)\b(phase|polar|exp|log10|sqrt|acos|asin|atan|cos|
|sin|tan|acosh|asinh|atanh|cosh|sinh|tanh|isfinite|isinf|isnan|log|rect|
|isclose)\(""", r"cmath.\1(", fn)

# so stuff like log(x)sin(x) works as expected
fn = re.sub(r"(\w|\))\s+(\w|\()", r"\1 * \2", fn)

# replace ^ with **
fn = re.sub(r"\^", r"**", fn)

fncode = compile(fn, "<string>", "eval")

c          =           input("c?              [0]:   ") or "0"
aspect     = max(float(input("aspect ratio?   [1]:   ") or 500), 0)
width      = max(  int(input("img width?      [500]: ") or 500), 1)
iterations = max(  int(input("iterations?     [32]:  ") or 32),  1)
crange     = max(float(input("c range?        [1.5]: ") or "1.5"), 0)
cellcount  = max(  int(input("cell count      [1]:   ") or 1), 1)
center     =           input("center?         [0 0]: ") or "0 0"
zoom       =     float(input("zoom?           [1]:   ") or "1")
colorscale = max(  int(input("gradient speed? [1]:   ") or 1), 1)
cutoff     = max(  int(input("escape radius?  [30]:  ") or 30), 0)

c = re.sub("\s", "", c)
c = re.sub("i", "j", c)
c = complex(c)

if crange <= 0:
    # ignored, but avoids divide by 0
    crange = 1
    cellcount = 1

height = int(width / aspect)
rowheight   = height / cellcount
columnwidth = width  / cellcount

cutoff = 30

inx = center.find(" ")
c_x = float(center[:inx])
c_y = float(center[(inx + 1):])
spread = 1 / zoom

graph = {
    "x": {
        "min": c_x - spread * aspect,
        "max": c_x + spread * aspect,
        "c": 0
    },
    "y": {
        "min": c_y - spread,
        "max": c_y + spread,
        "c": 0
    }
}

graph['x']['c'] = (graph['x']['max'] + graph['x']['min']) / 2
graph['y']['c'] = (graph['y']['max'] + graph['y']['min']) / 2

def write_pixel(r, g, b, f):
    f.write(bytes([r, g, b]))

# unix_time.ppm
fname_base = str(int(time.time()))
print("the time is " + fname_base)
with open("./info/" + fname_base + ".txt", "w") as out:
    out.write(
        "z_n+1 = " + orig_fn
        + "\nc = " + str(c.real) + " + " + str(c.imag) + "i"
        + "\nrendered area: ("
        + str(graph['x']['min']) + ", " + str(graph['y']['min']) + ") to ("
        + str(graph['x']['max']) + ", " + str(graph['y']['max']) + ")"
        + "\niterations: " + str(iterations)
        + "\n"
    )
    for row in range(1, cellcount + 1):
        localy = c.imag
        if cellcount != 1:
            localy += crange - 2 * crange * row / (cellcount - 1)

        for column in range(1, cellcount + 1):
            localc = c.real
            if cellcount != 1:
                localc += 2 * crange * column / (cellcount - 1) - crange
            out.write(str(column) + " " + str(row) + " = " + str(localc)
                + signstr(localy) + str(abs(localy)) + "i\n")

with open("./output/" + fname_base + ".ppm", "wb") as out:
    out.write(bytes("P6\n"
        + str(width)  + " " + str(height) + "\n255\n",
        encoding="ascii")
    )
    z: complex
    z_p: complex
    i: int
    color: float
    graphx: float
    row: int
    column: int
    localc: complex
    localy: float
    start = datetime.now()

    for y in range(0, height):

        # eta prediction, % done
        now = datetime.now()
        # Δt / % done = est. total time
        # ETT - Δt = est. time left
        doneamt = y / height
        if y != 0:
            eta = (now - start) * (1 / doneamt - 1)
            etastr = (
                "{: >6.3f}% done, eta ≈ {:02d}:{:02d}:{:02d}.{:03d}".format(
                100 * doneamt, # % complete
                math.floor(eta.seconds / 3600), # hours
                math.floor((eta.seconds % 3600) / 60), # minutes
                eta.seconds % 60, # seconds
                math.floor(eta.microseconds / 1000) # milliseconds
            ))
            print("{: <76}".format(etastr), end="\r")

        graphy = stgY(y % rowheight)
        row = math.floor(y / rowheight)
        localy = c.imag
        if cellcount != 1:
            localy += crange - 2 * crange * row / (cellcount - 1)

        for x in range(0, width):
            column = math.floor(x / columnwidth)

            localc = c.real + localy*1j

            if cellcount != 1:
                localc += 2 * crange * column / (cellcount - 1) - crange

            z_p = z = stgX(x % columnwidth) + graphy*1j
            i = 0
            # smooth coloring using exponents????
            color = math.exp(-abs(z))
            for i in range(0, iterations):
                z = eval_fn(z, localc)
                if cmath.isnan(z) or not cmath.isfinite(z):
                    # oh no
                    z = 1
                    break
                elif abs(z) > cutoff:
                    break
                z_p = z
                color += math.exp(-abs(z))

            color /= iterations
            write_pixel(
                int(255 * math.sin(color * colorscale * 9) ** 2),
                int(255 * math.sin(color * colorscale * 10) ** 2),
                int(255 * math.sin(color * colorscale * 11) ** 2),
                out
            )

    # print("", end="\n")
