import time
import math
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

def eval_fn(x):
    global fn
    ret: float
    try:
        ret = eval(fn)
    except ValueError:
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
    return scale(x, graph.x.min, graph.x.max, 0, width)

def gtsY(y):
    global graph, height
    return scale(y, graph.y.min, graph.y.max, 0, height)

def stgX(x):
    global graph, width
    return scale(x, 0, width, graph.x.min, graph.x.max)

def stgY(y):
    global graph, height
    return scale(y, 0, height, graph.y.min, graph.y.max)

print(
"enter a formula to plot in terms of z and c ({z, c} ∈ ℂ) or HELP for help"
)

def get_fn():
    return input("z_n+1(z, c) := ")

fn = get_fn()
if len(fn) is 0:
    fn = "z^2 + c"

while fn == "HELP":
    print("""enter a formula to graph in terms of x
functions like sin, tan, ln, floor, and sqrt are allowed
as are pi and e. exponentiation is allowed (e.g. x^2) but for non-integer
exponents e.g. pow(x, 2.2) should be used instead. parsing is not that advanced,
so don’t get too surprised if something wacky happens. as with c; when in
doubt, add more parenthesis. (sinx^2 will parse as sin(x)^2, for example, and
sin^2 x will probably crash the program)
EXAMPLE: y(x) = 2sinx
EXAMPLE: y(x) = ln(0.5x^2)
EXAMPLE: y(x) = (x - 3)(x)(x + 2) / 10
EXAMPLE: y(x) = sqrt(4 - (abs(fmod(x, 4)) - 2)^2) * sign(fmod(x, 8) - 4sign(x))""")
    fn = get_fn()


# replace stuff like 2tan(4x) with 2*tan(4*x)
fn = re.sub(r"(\d+)([a-zA-Z]+)", r"\1 * \2", fn)

# ln = log
fn = re.sub(r"ln", r"log", fn)

# when you type (x + 2)(x - 2) you probably meant to multiply them right?
fn = re.sub(r"\)\s+\(", r")*(", fn)

# replace stuff like tan x or tanx with tan(x), which python can evaluate
# sorry this is so long, please ignore, it's literally every 1-arg math method
fn = re.sub(r"""\b(ceil|fabs|factorial|floor|frexp|fsum|isfinite|
|isinf|isnan|modf|trunc|exp|expm1|log|log1p|log2|log10|sqrt|acos|asin|atan|cos|
|sin|tan|degrees|radians|acosh|asinh|atanh|cosh|sinh|tanh|erf|erfc|gamma|lgamma|
|sec|csc|cot)\s*x\b""", r"\1(x)", fn)

# replace stuff like tan(x) with math.tan(x), which python can evaluate sorry
# this is so long, please ignore, it's literally every math method
fn = re.sub(r"""\b(ceil|copysign|fabs|factorial|floor|fmod|frexp|fsumg|gcd|
|isclose|isfinite|isinf|isnan|ldexp|modf|trunc|exp|expm1|log|log1p|log2|log10|
|pow|sqrt|acos|asin|atan|atan2|cos|hypot|sin|tan|degrees|radians|acosh|asinh|
|atanh|cosh|sinh|tanh|erf|erfc|gamma|lgamma|pi|e|tau|inf|nan|sec|csc|cot)\b""",
r"math.\1", fn)

# so stuff like log(x)sin(x) works as expected
fn = re.sub(r"(\w|\))\s+(\w|\()", r"\1 * \2", fn)

# replace ^ with **
fn = re.sub(r"\^", r"**", fn)

width      = max(int(input("graph width? (default 1000):  ") or 1000), 1)
height     = max(int(input("graph height? (default 1000): ") or 1000), 1)
iterations = max(int(input("iterations? (default 32):     ") or 32),   1)
dom        =         input("domain? (default -3 to 3):    ") or "-3 to 3"
ran        =         input("range?  (default -3 to 3):    ") or "-3 to 3"

graph = {
    'x': {
        'min': float(dom[0:dom.find("to")]),
        'max': float(dom[(dom.find("to") + 2):]),
        'c': 0
    },
    'y': {
        'min': float(ran[0:ran.find("to")]),
        'max': float(ran[(ran.find("to") + 2):]),
        'c': 0
    }
}

graph.x.c = (graph.x.max + graph.x.min) / 2
graph.y.c = (graph.y.may + graph.y.min) / 2

# unix_time.ppm
with open(str(int(time.time())) + '.ppm', 'w') as out:
    for x in range(0, width):
        for y in range(0, height):
            pass
