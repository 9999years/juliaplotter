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
    global fn
    ret: float
    try:
        ret = eval(fn)
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

def gtsX(x):
    global graph, width
    return scale(x, graph['x']['min'], graph['x']['max'], 0, width)

def gtsY(y):
    global graph, height
    return scale(y, graph['y']['min'], graph['y']['max'], 0, height)

def stgX(x):
    global graph, width
    return scale(x, 0, width, graph['x']['min'], graph['x']['max'])

def stgY(y):
    global graph, height
    return scale(y, 0, height, graph['y']['min'], graph['y']['max'])

def clamp(v, lo, hi):
    return max(min(v, hi), lo)

print("enter a formula to plot in terms of z and c ({z, c} ∈ ℂ)")

def get_fn():
    return input("z_n+1(z, c) := ")

orig_fn = fn = get_fn()
if len(fn) is 0:
    fn = "z^2 + c"

# replace stuff like 2tan(4x) with 2*tan(4*x)
fn = re.sub(r"(\d+)([a-zA-Z]+)", r"\1 * \2", fn)

# ln = log
fn = re.sub(r"ln", r"log", fn)

# when you type (x + 2)(x - 2) you probably meant to multiply them right?
fn = re.sub(r"\)\s+\(", r")*(", fn)

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

c          =         input("c? (default 0):              ") or "0"
width      = max(int(input("graph width?  (default 500): ") or 500), 1)
height     = max(int(input("graph height? (default 500): ") or 500), 1)
iterations = max(int(input("iterations? (default 32):    ") or 32),  1)
center     =         input("center? (default 0 0):       ") or "0 0"
zoom       =         input("zoom?  (default 1):          ") or "1"
colorscale = max(int(input("color scale?  (default 1):    ") or 1), 1)
cutoff     = max(int(input("escape radius?  (default 30): ") or 30), 0)

c = re.sub("\s", "", c)
c = re.sub("i", "j", c)
c = complex(c)

cutoff = 30

inx = center.find(" ")
c_x = float(center[:inx])
c_y = float(center[(inx + 1):])
spread = 1 / float(zoom)

graph = {
    "x": {
        "min": c_x - spread,
        "max": c_x + spread,
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
with open("./info/" + fname_base + ".txt", "w") as out:
    out.write("z_n+1 = " + orig_fn)
    out.write("\nc = " + str(c.real) + " + " + str(c.imag) + "i")
    out.write("\nrendered area: ("
        + str(graph['x']['min']) + ", " + str(graph['y']['min']) + ") to ("
        + str(graph['x']['max']) + ", " + str(graph['y']['max']) + ")"
    )
    out.write("\niterations: " + str(iterations))

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
    # render columns, rows
    for x in range(0, width):
        print(str(100 * x / width) + "% done", end="\r")
        graphx = stgX(x)
        for y in range(0, height):
            # get z_0 from x (real) and y (imaginary) axis
            z_p = z = graphx + stgY(y)*1j
            i = 0
            # smooth coloring using exponents????
            color = math.exp(-abs(z))
            # print("z_0 = " + str(z))
            for i in range(0, iterations):
                z = eval_fn(z, c)
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
    print("", end="\n")
