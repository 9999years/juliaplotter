# for testing shader algs

import time
import math
import codecs
import sys

# linear map val∈[valmin, valmax] |→ out∈[outmin, outmax]
def scale(val,  valmin,  valmax,  outmin,  outmax):
    return (
        (val - valmin) / (valmax - valmin)
        * (outmax - outmin) + outmin
    )

def write_pixel(r, g, b, f):
    f.write(bytes([r, g, b]))

width = 500
height = 500

# multiples of π:
# 3.1415926535
# 6.283185307
# 9.4247779605
# 12.566370614
# 15.7079632675
# 18.849555921
# multiples of π/2:
# 1.57079632675
# 3.1415926535
# 4.71238898025
# 6.283185307
# 7.85398163375
# 9.4247779605
# 10.99557428725
# 12.566370614
# 14.13716694075
# 15.7079632675
# 17.27875959425
# 18.849555921

# unix_time.ppm
fname_base = str(int(time.time()))
print(fname_base)
with open("./output/test" + fname_base + ".ppm", "wb") as out:
    out.write(bytes("P6\n"
        + str(width)  + " " + str(height) + "\n255\n",
        encoding="ascii")
    )
    color: float
    r: int
    g: int
    b: int
    # render columns, rows
    for x in range(0, width):
        color = 1.2 * x / width
        r = int(255 * math.sin(color * 9) ** 2)
        g = int(255 * math.sin(color * 10) ** 2)
        b = int(255 * math.sin(color * 11) ** 2)
        for y in range(0, height):
            write_pixel(r, g, b, out)
