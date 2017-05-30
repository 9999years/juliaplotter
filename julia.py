# for eta
from datetime import datetime
# for unix timestamp filenames
import time
# floor mostly, some assorted stuff
import math
# misc for actually evaluating the julia equations
import cmath
# turning user input into an eval-able formula
import re
# file writing
import codecs
# args? i actually might not need this
import argparse

# why arent these in math in the first place?
# we need them to make user formulae sensible
# maybe i should make a pull request...
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

# evaluate user input formula, compiled to bytecode
def eval_fn(z, c):
    global fncode
    try:
        return eval(fncode)
    except (ValueError, OverflowError, ZeroDivisionError):
        # negative number in a log or root probably
        # probably a user error
        return float('nan')

# linear map val∈[valmin, valmax] |→ out∈[outmin, outmax]
def scale(val, valmin, valmax, outmin, outmax):
    return (
        (val - valmin) / (valmax - valmin)
        * (outmax - outmin) + outmin
    )

# linear map val∈[0, valmax] |→ out∈[outmin, outmax]
def zeroscale(val, valmax, outmin, outmax):
    return (val / valmax) * (outmax - outmin) + outmin

# screen coords to graph coords
def stgX(x):
    global graph, colwidth
    return zeroscale(x % colwidth, colwidth,
            graph['x']['min'], graph['x']['max'])

def stgY(y):
    global graph, rowheight
    return zeroscale(y % rowheight, rowheight,
        graph['y']['max'], graph['y']['min'])

def clamp(v, lo, hi):
    return max(min(v, hi), lo)

# printing complex numbers
def signstr(num):
    if math.sign(num) == -1:
        return '-'
    else:
        return '+'

def strcomplex(num):
    return '{:8g} {} {:<8g}i'.format(
        num.real,
        signstr(num.imag),
        abs(num.imag))

def process_fn(fn):
    # replace stuff like 2tan(4x) with 2*tan(4*x)
    fn = re.sub(r'(\d+)([a-zA-Z]+)', r'\1 * \2', fn)

    # ln = log
    fn = re.sub(r'ln', r'log', fn)

    # when you type (x + 2)(x - 2) you probably meant to multiply them right?
    fn = re.sub(r'\)\s*\(', r')*(', fn)

    fn = re.sub(r'π', r'pi', fn)

    fn = re.sub(r'(?<!cmath\.)\b(pi|e|tau|inf|infj|nan|nanj)\b',
        r'cmath.\1', fn)

    # (?<! ...)

    # sinz, sin z, cos c, etc.
    fn = re.sub(r'''(?<!cmath\.)\b(phase|polar|exp|log10|sqrt|acos|asin|atan|
    |cos|sin|tan|acosh|asinh|atanh|cosh|sinh|tanh|isfinite|isinf|isnan|log|
    |rect)\s*([zc])\b''', r'\1(\2)', fn)

    # sin(z) ...
    fn = re.sub(r'''(?<!cmath\.)\b(phase|polar|exp|log10|sqrt|acos|asin|atan|
    |cos|sin|tan|acosh|asinh|atanh|cosh|sinh|tanh|isfinite|isinf|isnan|log|
    |rect|isclose)\(''', r'cmath.\1(', fn)

    # so stuff like log(x)sin(x) works as expected
    fn = re.sub(r'(\w|\))\s+(\w|\()', r'\1 * \2', fn)

    # replace ^ with **
    fn = re.sub(r'\^', r'**', fn)

    return fn

def desc(description):
    if description[0] == '\n':
        return description[1:].replace('\n', ' ')
    else:
        return description.replace('\n', ' ')

parser = argparse.ArgumentParser(
    description='Render an arbitrary Julia set or a grid of Julia sets'
    'with differing constants c from a user-entered equation.'
)

parser.add_argument('--fn', '-f', metavar='zₙ₊₁', type=str,
    default='z**2 + c', help='''The Julia set's function for iteration.''')

parser.add_argument('-c', metavar='constant', type=str,
    default='0', help=desc('''
The constant c for the function zₙ₊₁(z, c). Enter `random` to select a random
value for c.'''))

parser.add_argument('-a', '--aspect', metavar='aspect', type=float,
    default=1.0, help=desc('''
The output image's w/h aspect ratio.  Ex.: -a 2 implies an image twice as wide
as it is tall.'''))

parser.add_argument('-w', '--width', metavar='width', type=int,
    default='500', help='''The output image\'s width.''')

parser.add_argument('-i', '--iterations', metavar='iterations', type=int,
    default=32, help='The iterations to calculate the set to.')

parser.add_argument('-r', '--c-range', metavar='c-range', type=float,
    default=1.5, help=desc('''
The range of c values to use --- only relevant if the cell count option is used
to render a grid of sets; the c values for each sets will range from (c_r -
crange, c_i - crange·i) to (c_r + crange, c_i + crange·i), where c_r and c_i
are the real and imaginary components of the constant supplied with -c.'''))

parser.add_argument('-e', '--center', metavar='center', type=float,
    default=[0, 0], nargs=2, help=desc('''
The coordinate the graph is centered around, entered as two floats separated by
a space. (Not a comma! No parenthesis!)'''))

parser.add_argument('-n', '--cell-count', metavar='cell count', type=int,
    default=1, help=desc('''
The number of rows and columns to render. A cell count of 1 will render a
single set, and other values will render grids of Julia sets. The different
values of c are determined by --c-range or -r.'''))

parser.add_argument('-z', '--zoom', metavar='zoom', type=float,
    default=1, help=desc('''
How zoomed in the render is. The distance between the center-point and the top
/ bottom of the rendered area is 1 / zoom. Larger values of will produce a more
zoomed-in image, smaller values (<1) will produce a more zoomed-out image.'''))

parser.add_argument('-g', '--gradient', metavar='gradient speed', type=float,
    default=1, help=desc('''
The plotter colors images by smoothly interpolating the orbit escape times for
each value of z₀ in the, governed by a sine function. This option adds a
multiplier within the sine function to increase the oscillation speed, which
may help to enhance details in lightly colored images.'''))

parser.add_argument('-u', '--cutoff', '--escape-radius',
    metavar='escape', type=float, default=30, help=desc('''
The orbit escape radius --- how large |zₙ| must be before it's considered to
have diverged. Usually ≈ 30 for Julia sets, 2 for the Mandelbrot set.'''))

parser.add_argument('-o', '--output', metavar='directory', type=str,
    default='./output/', help='Output directory to write images to.')

parser.add_argument('--info-dir-name',
    metavar='directory', type=str, default='info', help=desc('''
Directory to write information files to. Is always a first-level directory
within the output directory, and changing it will probably mess up HTML
output.'''))

parser.add_argument('--no-info', action='store_false',
    help='''Don't write the HTML info file.''')

parser.add_argument('--no-render', action='store_true', help=desc('''
Generates appropriate HTML information files but doesn't render an
image (useful with `--filename` if the info for an image has been lost to
the sands of time...)'''))

parser.add_argument('--no-progress', action='store_false', help=desc('''
Don't output progress percentage and finish ETA. May increase
performance.'''))

parser.add_argument('--filename', metavar='pathspec', type=str, help=desc('''
Filename base for the output image. Defaults to the Unix timestamp. Relative to
the output directory. Shouldn't include extensions.'''))

parser.add_argument('--no-convert', action='store_false', help=desc('''
Don't shell out to `magick` to convert the .ppm to a .png after rendering.'''))

parser.add_argument('--no-open', action='store_false', help=desc('''
Don't open HTML output in a browser after completing rendering.'''))

args = parser.parse_args()

aspect     = args.aspect
c          = args.c
crange     = args.c_range
cellcount  = args.cell_count
center     = args.center
cutoff     = args.cutoff
fn         = args.fn
iterations = args.iterations
colorscale = args.gradient * iterations / 32
width      = args.width
zoom       = args.zoom
info_dir   = args.info_dir_name
out_dir    = args.output
no_render  = args.no_render
fname = args.filename
show_prog  = args.no_progress
write_info = args.no_info
convert    = args.no_convert
open_html  = args.no_open

def zero_warning(var, extra=''):
    print('WARNING: ' + var + ' is ZERO. Ignoring. ' + extra)

def neg_warning(var, extra=''):
    print('WARNING: ' + var + ' is NEGATIVE.'
        + 'This makes no sense, using absolute value instead. ' + extra)

# validation
if aspect == 0:
    zero_warning('aspect')
    aspect = 1
elif aspect < 0:
    neg_warning('aspect')
    aspect = abs(aspect)

if crange == 0:
    zero_warning('c range', extra='Rendering only one cell.')
    crange = 1
    cellcount = 1
elif crange < 0:
    neg_warning('c range')
    crange = abs(crange)

if cellcount == 0:
    zero_warning('cell count')
    cellcount = 1
elif cellcount < 0:
    neg_warning('cell count')
    cellcount = abs(cellcount)

if cutoff == 0:
    zero_warning('escape radius')
    cutoff = 30
elif cutoff < 0:
    neg_warning('escape radius')
    cutoff = abs(cutoff)

if iterations == 0:
    zero_warning('iterations')
    iterations = 1
elif iterations < 0:
    neg_warning('iterations')
    iterations = abs(iterations)

if colorscale == 0:
    zero_warning('gradient scale')
    colorscale = iterations / 32
elif colorscale < 0:
    neg_warning('gradient scale')
    colorscale = abs(colorscale)

if zoom == 0:
    zero_warning('zoom')
    zoom = 1
elif zoom < 0:
    neg_warning('zoom')
    zoom = abs(zoom)

orig_fn = fn
fn = process_fn(fn)

fncode = compile(fn, '<string>', 'eval')

if c.lower() == 'random':
    from random import random
    c = 2 * random() - 1 + 2j * random() - 1j
    print('c = {}'.format(strcomplex(c)))
else:
    c = re.sub('\s', '', c)
    c = re.sub('i', 'j', c)
    c = complex(c)

height = int(width / aspect)
rowheight = height / cellcount
colwidth  = width  / cellcount

cutoff = 30

c_x = center[0]
c_y = center[1]
spread = 1 / zoom

graph = {
    'x': {
        'min': c_x - spread * aspect,
        'max': c_x + spread * aspect,
        'c': 0
    },
    'y': {
        'min': c_y - spread,
        'max': c_y + spread,
        'c': 0
    }
}

graph['x']['c'] = (graph['x']['max'] + graph['x']['min']) / 2
graph['y']['c'] = (graph['y']['max'] + graph['y']['min']) / 2

cgrid = [[0 for y in range(cellcount)] for x in range(cellcount)]
#cgrid[x][y]
yticks = [int((y + 1) * rowheight) for y in range(cellcount - 1)]
xticks = [int((x + 1) * colwidth ) for x in range(cellcount - 1)]

for row in range(cellcount):
    localy = c.imag
    if cellcount != 1:
        localy += crange - 2 * crange * row / (cellcount - 1)

    for col in range(cellcount):
        localc = c.real
        if cellcount != 1:
            localc += 2 * crange * col / (cellcount - 1) - crange
        cgrid[col][row] = localc + localy*1j

# unix_time.ppm
fname = fname or str(int(time.time()))
print('the time is ' + fname)

if write_info:
    with open('./starttemplate.html') as template_start, \
        open('./endtemplate.html') as template_end, \
        open(out_dir + '/' + info_dir + '/' + fname + '.html',
            'w', encoding='utf-8') as out:
        targets_str = ''
        out.write(template_start.read())
        if cellcount > 1:
            out.write('<map name="juliagrid" id="juliamap">')
            for row in range(cellcount):
                for col in range(cellcount):
                    out.write(
                        '<area shape="rect" coords="' + ','.join(
                        [str(int(colwidth  *  col)) #top left
                        ,str(int(rowheight *  row))
                        ,str(int(colwidth  * (col + 1))) #bottom right
                        ,str(int(rowheight * (row + 1)))])
                        + f'" href="#{col + 1}-{row + 1}">\n'
                    )
                    targets_str += (
                        f'<div id="{col + 1}-{row + 1}">column {col + 1}, '
                        f'row {row + 1}: c = {strcomplex(cgrid[col][row])}'
                        '<p><code>python julia.py '
                        f'-f "{orig_fn}" -c "{strcomplex(cgrid[col][row])}" '
                        f'-i {iterations} -w {width} -a {aspect} '
                        f'-e {args.center[0]} {args.center[1]} -z {zoom} '
                        f'-g {colorscale} -u {cutoff}</code></div>\n'
                    )
            out.write('</map>\n')
        else:
            targets_str = (
                '<p><code>python julia.py '
                f'-f "{orig_fn}" -c "{strcomplex(cgrid[0][0])}" '
                f'-i {iterations} -w {width} -a {aspect} '
                f'-e {args.center[0]} {args.center[1]} -z {zoom} '
                f'-g {colorscale} -u {cutoff}</code>\n'
            )

        def tr(one, two):
            return f'<tr><td>{one}</td><td>{two}</td></tr>'

        out.write(
            '<img ' + ('usemap="#juliagrid" ' if cellcount > 1 else '')
            + f'src="../{fname}.png" id="julia">\n' +
            ('<div class="targets">' if cellcount > 1 else '')
            + targets_str +
            (('</div><p>click on a cell to see the constant and the '
            'command-line invocation used to render it!') if cellcount > 1
            else '')
            + '<table class="render-info">'
            + tr('zₙ₊₁', orig_fn)
            + tr('c', strcomplex(c))
            + tr('rendered area',
                f'({graph["x"]["min"]}, {graph["y"]["min"]}i)'
                f' to ({graph["x"]["max"]}, {graph["y"]["max"]}i)')
            + tr('center', f'({graph["x"]["c"]}, {graph["y"]["c"]}i)')
            + tr('zoom', f'{zoom}×')
            + tr('gradient speed', colorscale)
            + tr('escape radius', cutoff)
            + tr('iterations', iterations)
            + tr('c range', f'{crange:g}')
            + '</table>'
            + template_end.read())

if no_render:
    exit()

def write_pixel(r, g, b, f):
    f.write(bytes([r, g, b]))

def strtimedelta(time):
    return (f'{math.floor(time.seconds / 3600):02d}:'
        f'{math.floor((time.seconds % 3600) / 60):02d}:'
        f'{time.seconds % 60:02d}.{math.floor(time.microseconds / 1000):03d}')

with open(out_dir + '/' + fname + '.ppm', 'wb') as out:
    out.write(bytes('P6\n{} {}\n255\n'.format(width, height),
        encoding='ascii'))
    z: complex
    z_p: complex
    i: int
    color: float
    graphy: float
    row: int = 0
    col: int = 0
    start = datetime.now()

    for y in range(0, height):

        if show_prog:
            # eta prediction, % done
            now = datetime.now()
            # Δt / % done = est. total time
            # ETT - Δt = est. time left
            doneamt = y / height
            if y != 0:
                eta = (now - start) * (1 / doneamt - 1)
                print('{: <76}'.format(
                    f'{100 * doneamt: >6.3f}% done, eta ≈ {strtimedelta(eta)}'
                ), end='\r')

        graphy = stgY(y) * 1j
        if row != cellcount - 1 and y == yticks[row]:
            row += 1

        col = 0
        for x in range(0, width):
            if col != cellcount - 1 and x == xticks[col]:
                col += 1

            z_p = z = stgX(x) + graphy
            i = 0
            # smooth coloring using exponents????
            color = math.exp(-abs(z))
            for i in range(0, iterations):
                z = eval_fn(z, cgrid[col][row])
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
                int(255 * math.sin(color * colorscale * 9 ) ** 2),
                int(255 * math.sin(color * colorscale * 10) ** 2),
                int(255 * math.sin(color * colorscale * 11) ** 2),
                out
            )

    end = datetime.now()
    print(f'Done! Completed in {strtimedelta(end - start)}')

if convert:
    # convert ppm to png (if requested, by default on)
    from subprocess import run
    import os
    print(f'Converting {fname}.ppm to {fname}.png')
    run(['magick', 'mogrify', '-format', 'png', f'{out_dir}/{fname}.ppm'])
    os.remove(f'{out_dir}/{fname}.ppm')
if open_html:
    import webbrowser
    from urllib.request import pathname2url
    webbrowser.open('file:'
        + pathname2url(os.path.abspath(
            f'{out_dir}/{info_dir}/{fname}.html'
        ))
    )
