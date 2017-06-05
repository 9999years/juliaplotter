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
# warnings
import sys
# making directories, deleting the .ppm after conversion, resolving paths
import os

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

# https://stackoverflow.com/a/14981125/5719760
# sys module for stderr
# import sys
# def error(*args, **kwargs):
    # print('ERROR: ', *args, file=sys.stderr, **kwargs)

def warn(*args, **kwargs):
    print('WARNING: ', *args, file=sys.stderr, **kwargs)

# evaluate user input formula, compiled to bytecode
def eval_fn(z, c):
    global fncode
    try:
        return eval(fncode)
    except (ArithmeticError, ValueError, OverflowError, ZeroDivisionError):
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
    if math.sign(num.real) == -1:
        return '-'
    else:
        return '+'

#something like 0.0 - 3.2i
def strcomplex(num):
    return f'{num.real:8g} {signstr(num.imag)} {abs(num.imag):<8g}i'.strip()

def process_fn(fn):
    global orig_fn
    if fn.lower() == 'random':
        # general case:
        # https://gist.github.com/9999years/7d9430d08cba96c928b5631293dde33e
        import random

        def get_coef(imag=True):
            max_coef = 20
            coef = (random.random() - 0.5) * max_coef
            if imag:
                coef += get_coef(imag=False) * 1j
            return coef

        def term(deg, variables, consts, imag):
            deg = random.randint(0, deg)
            ret = '('
            var = 'z'
            coef = get_coef(imag)
            if deg != 0:
                ret += f'{var}'
                if deg > 1:
                    ret += f'^{deg} + '
                else:
                    ret += ' + '
            ret += f'({coef:.3f}){random.choice(consts)})'
            return (deg, ret)

        def poly(deg, variables, consts, imag):
            ret = ""
            i = 0
            while i < random.randrange(1, deg):
                (degtmp, rtmp) = term(deg, variables, consts, imag)
                i += degtmp
                ret += rtmp
            return ret

        def rational(deg, variables, consts, imag):
            ret = poly(deg, variables, consts, imag)
            if random.random() > 0.3:
                ret += f'/({poly(deg * 2, variables, consts, imag)})'
            return ret

        degree    = 3
        variables = ['z']
        consts    = ['c', '']
        imag      = True

        fn = orig_fn = rational(degree, variables, consts, imag)
        orig_fn = re.sub(r'j', r'i', fn)
        print(f'f(z, c) = {orig_fn}')
        fn = re.sub(r'\^', r'**', fn)
        fn = re.sub(r'([zc)])([zc(])', r'\1 * \2', fn)
        return fn

    # imaginary numbers are j
    fn = re.sub(r'(\d|\b)i\b', r'\1j', fn)

    # ln = log
    fn = re.sub('ln', 'log', fn)

    # when you type (x + 2)(x - 2) you probably meant to multiply them right?
    fn = re.sub(r'\)\s*\(', r')*(', fn)

    fn = re.sub('π', 'pi', fn)

    fn = re.sub(r'(?<!cmath\.)\b(pi|e|tau|inf|infj|nan|nanj)\b',
        r'cmath.\1', fn)

    # sinz, sin z, cos c, etc.
    fn = re.sub(r'''(?<!cmath\.)(?:(?<=\d)|(?<=\b))(phase|polar|exp|log10|sqrt|
    |acos|asin|atan| |cos|sin|tan|acosh|asinh|atanh|cosh|sinh|tanh|isfinite|
    |isinf|isnan|log|rect)\s*([zc])\b''', r'cmath.\1(\2)', fn)

    # sin(z) ...
    fn = re.sub(r'''(?<!cmath\.)(?:(?<=\d)|(?<=\b))(phase|polar|exp|log10|sqrt|acos|asin|atan|
    |cos|sin|tan|acosh|asinh|atanh|cosh|sinh|tanh|isfinite|isinf|isnan|log|
    |rect|isclose)\(''', r'cmath.\1(', fn)

    # replace stuff like 2tan(4x) with 2*tan(4*x)
    # (?=j\b) excludes imaginary numbers
    fn = re.sub(r'(\d+)(?!j\b)\s*([a-zA-Z]+)', r'\1 * \2', fn)

    # 3 z c, 2.5c z
    fn = re.sub(r'([zc])\s+([zc])', r'\1 * \2', fn)

    # so stuff like x(x - 3) works as expected
    fn = re.sub(r'([zc])\s+\(', r'\1 * (', fn)
    # so stuff like log(x)sin(x) works as expected
    fn = re.sub(r'\)\s*(\w)', r') * \1', fn)

    # replace ^ with **
    fn = re.sub(r'\^', r'**', fn)

    z = 0
    c = 0
    try:
        eval(fn)
    except (ArithmeticError, ValueError):
        warn('Evaluation of function fails at z = 0, c = 0! This might be a '
            'symptom of a larger problem, or simply a harmless asymptote. '
            'Continuing execution.')
    except NameError:
        warn('Uh-oh, something is pretty seriously wrong in the function you '
            'gave me. Continuing execution, but this is *probably* going to '
            'explode in a moment.')
        print(f'Processed function: {fn}')

    return fn

# so i can nicely format my help messages and not worry about weird
# formatting for the user (sorry!)
def desc(description):
    if description[0] == '\n':
        return description[1:].replace('\n', ' ')
    else:
        return description.replace('\n', ' ')

# set up arguments
parser = argparse.ArgumentParser(
    description='Render an arbitrary Julia set or a grid of Julia sets'
    'with differing constants c from a user-entered equation.'
)

parser.add_argument('--fn', '-f', metavar='zₙ₊₁', type=str,
    default=None, help=desc('''
The Julia set's function for iteration. Enter `random` to generate a random
complex rational function P(z)/Q(z), where P(z) and Q(z) are complex polynomials
of maximum degree 3 and 6, respectively.'''))

parser.add_argument('-c', '--constant', metavar='constant', type=None,
    default='0', help=desc('''
The constant c for the function zₙ₊₁(z, c). Enter `random` to select a random
value for c. Default: 0 + 0i'''))

parser.add_argument('-a', '--aspect', metavar='aspect', type=float,
    default=1.0, help=desc('''
The output image's w/h aspect ratio.  Ex.: -a 2 implies an image twice as wide
as it is tall. Default: 1.0'''))

parser.add_argument('-w', '--width', metavar='width', type=int,
    default='500', help='''The output image\'s width.''')

parser.add_argument('-i', '--iterations', metavar='iterations', type=int,
    default=32, help='The iterations to calculate the set to.')

parser.add_argument('-r', '--c-range', metavar='c-range', type=float,
    default=1.5, help=desc('''
The range of c values to use --- only relevant if the cell count option is used
to render a grid of sets; the c values for each sets will range from (c_r -
crange, c_i - crange·i) to (c_r + crange, c_i + crange·i), where c_r and c_i
are the real and imaginary components of the constant supplied with -c. Default:
1.5'''))

parser.add_argument('-n', '--cell-count', metavar='cell count', type=int,
    default=1, help=desc('''
The number of rows and columns to render. A cell count of 1 will render a
single set, and other values will render grids of Julia sets. The different
values of c are determined by --c-range or -r. Default: 1'''))

parser.add_argument('-e', '--center', metavar='center', type=float,
    default=[0, 0], nargs=2, help=desc('''
The coordinate the graph is centered around, entered as two floats separated by
a space. (Not a comma! No parenthesis! It's technically two separate arguments
consumed by one option.) Default: 0 0'''))

parser.add_argument('-z', '--zoom', metavar='zoom', type=float,
    default=1, help=desc('''
How zoomed in the render is. The distance between the center-point and the top
/ bottom of the rendered area is 1 / zoom. Larger values of will produce a more
zoomed-in image, smaller values (<1) will produce a more zoomed-out image.
Default: 1'''))

parser.add_argument('-g', '--gradient', metavar='gradient speed', type=float,
    default=1, help=desc('''
The plotter colors images by smoothly interpolating the orbit escape times for
each value of z₀ in the, governed by a sine function. This option adds a
multiplier within the sine function to increase the oscillation speed, which
may help to enhance details in lightly colored images. Default: 1.0'''))

parser.add_argument('-u', '--cutoff', '--escape-radius',
    metavar='escape', type=float, default=30, help=desc('''
The orbit escape radius --- how large |zₙ| must be before it's considered to
have diverged. Usually ≈ 30 for Julia sets, 2 for the Mandelbrot set. Default:
30.0'''))

parser.add_argument('-o', '--output', metavar='directory', type=str,
        default='./output/', help=desc('''
Output directory to write images to. Default: ./output/'''))

parser.add_argument('--info-dir',
    metavar='directory', type=str, default='./info/', help=desc('''
Directory to write information files to, relative to the output directory. If
it’s not a first-level directory within the output directory, HTML output will
look funny. Default: ./info/'''))

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
Filename base for the output image. Relative to the output directory. Shouldn't
include extensions. Defaults: The current Unix timestamp'''))

parser.add_argument('--no-convert', action='store_false', help=desc('''
Don't shell out to `magick` to convert the .ppm to a .png after rendering.'''))

parser.add_argument('--no-open', action='store_false', help=desc('''
Don't open HTML output in a browser after completing rendering.'''))

parser.add_argument('-s', '--silent', action='store_true', help=desc('''
Don't log info, show progress, convert the .ppm to a .png, or open the file when
finished rendering.  Equivalent to `--no-open --no-convert --no-progress
--no-info`.'''))

# parse arguments, extract variables
args = parser.parse_args()

aspect     = args.aspect
c          = args.constant
crange     = args.c_range
cellcount  = args.cell_count
center     = args.center
cutoff     = args.cutoff
fn         = args.fn
iterations = args.iterations
# correct coloring for iterations so that coloring is same at different
# iterations. stumbled upon this fix by accident (guessing), no clue why it
# works honestly
colorscale = args.gradient * iterations / 32
width      = args.width
zoom       = args.zoom
info_dir   = args.info_dir
out_dir    = args.output
no_render  = args.no_render
fname      = args.filename
show_prog  = args.no_progress
write_info = args.no_info
convert    = args.no_convert
open_html  = args.no_open

if args.silent:
    show_prog = write_info = convert = open_html = False

def zero_warning(var, extra=''):
    warn(var + ' is ZERO. Ignoring. ' + extra)

def neg_warning(var, extra=''):
    warn(var + ' is NEGATIVE.'
        + 'This makes no sense, using absolute value instead. ' + extra)

# validate arguments to prevent weird nonsense / errors
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

if c is None and fn is None:
    fn = 'z^2 + c'
    c = 'random'
elif c is None:
    c = '0 + 0i'
elif fn is None:
    fn = 'z^2 + c'

# save original function, process to be usable and compile
orig_fn = fn
if orig_fn.lower() == 'random':
    c = 'random'
fn = process_fn(fn)

try:
    fncode = compile(fn, '<string>', 'eval')
except (SyntaxError, ValueError):
    raise Exception('Invalid function. This might be my fault or yours.\n'
        f'Processed equation: {fn}')
    exit()

# if the user wants a random c, generate one and tell them about it
if c.lower() == 'random':
    from random import random
    c = 2 * random() - 1 + 2j * random() - 1j
    print('c = {}'.format(strcomplex(c)))
else:
    # otherwise, parse the user's c and store in `c`
    c = re.sub('\s', '', c)
    c = re.sub('i', 'j', c)
    c = complex(c)

# aspect = width / height ⇒ height = width / aspect
# should be an int so the iteration stuff doesn't explode
height = int(width / aspect)
# more accuracy if we don't int-ify these (probably)
rowheight = height / cellcount
colwidth  = width  / cellcount

# two args for center variables
c_x = center[0]
c_y = center[1]

# so that a higher zoom = a smaller graph
spread = 1 / zoom

graph = {
    'x': {
        # make sure the graph isn't stretched
        'min': c_x - spread * aspect,
        'max': c_x + spread * aspect,
        'c':   c_x
    },
    'y': {
        'min': c_y - spread,
        'max': c_y + spread,
        'c':   c_y
    }
}

# it seems terrible that there's no better way to initialize a 2d array!
# anyways this reads backwards and keys are addressed as cgrid[x][y]
if cellcount is 1:
    # avoids a divide by zero
    cgrid = [[c]]
else:
    # tick values range from c - crange to c + crange
    # y axis is inverted so imag and real components arent the same
    cgrid = [[
           c.real - crange + 2 * crange * col / (cellcount - 1)
        + (c.imag + crange - 2 * crange * row / (cellcount - 1))*1j
        for row in range(cellcount)]
        for col in range(cellcount)]

# when we should step the row / col counters
yticks = [int((y + 1) * rowheight) for y in range(cellcount - 1)]
xticks = [int((x + 1) * colwidth ) for x in range(cellcount - 1)]

# unix timestamp is the base filename
fname = fname or str(int(time.time()))
# output it for future reference (if needed)
print('the time is ' + fname)

# make sure the directories we want to write to exist
# i know i said the info dir has to be in the output dir but if you
# really wanna be mean and fuck it up that's fine (the css and js links are
# gonna break though unless you modify starttemplate.html and endtemplate.html)
if not os.path.exists(f'./{out_dir}/'):
    os.makedirs(f'./{out_dir}/')
if not os.path.exists(f'./{out_dir}/{info_dir}/'):
    os.makedirs(f'./{out_dir}/{info_dir}/')

# write html if requested
# perhaps interestingly, we do this before actually generating the image
if write_info:
    # we're going to need to output render commands a lot --- abstract it
    # we need to vary the given c and toggle showing the cell count, so
    # introduce options for that
    def generateCLIinvocation(c, showcells=False):
        global orig_fn, iterations, width, aspect, args, \
            zoom, colorscale, cutoff, cellcount
        return ('./julia.py '
            f'-f "{orig_fn}" -c "{strcomplex(c)}" '
            f'-i {iterations} -w {width} -a {aspect} '
            f'-e {args.center[0]} {args.center[1]} -z {zoom} '
            f'-g {colorscale} -u {cutoff}'
            + (f' -r {crange:g} -n {cellcount}' if showcells else ''))

    # get template information. js and css are just linked to, although they
    # rely on writing the html two directories below this directory (usually
    # in ./output/info)
    with open('./starttemplate.html', encoding='utf-8') as template_start, \
        open('./endtemplate.html', encoding='utf-8') as template_end, \
        open(out_dir + '/' + info_dir + '/' + fname + '.html',
            'w', encoding='utf-8') as out:
        # targets is a string containing info to replicate the render or
        # parts of the render. it contains <div>s for each row / col wrapped
        # in a div.targets if cellcount > 1 with the cli invocation, or just
        # one not in a div.targets if cellcount = 1
        targets = ''
        end_script = ''
        out.write(template_start.read())
        if cellcount > 1:
            out.write('<map name="juliagrid" id="juliamap">')
            targets += '<div class="targets">'
            for row in range(cellcount):
                for col in range(cellcount):
                    # write an image map with an area for each row/col. then
                    # when clicked, the relevant info from targets is shown
                    # (displayed in a :target css selector, hidden by default)
                    out.write(
                        '<area shape="rect" coords="' + ','.join(
                        [str(int(colwidth  *  col)) #top left
                        ,str(int(rowheight *  row))
                        ,str(int(colwidth  * (col + 1))) #bottom right
                        ,str(int(rowheight * (row + 1)))])
                        + f'" href="#{col + 1}-{row + 1}">\n'
                    )
                    targets += (
                        f'<div id="{col + 1}-{row + 1}">column {col + 1}, '
                        f'row {row + 1}: c = {strcomplex(cgrid[col][row])}'
                        '<p><code>'
                        + generateCLIinvocation(cgrid[col][row])
                        + '</code></div>\n'
                    )
            out.write('</map>\n<img usemap="#juliagrid" ')
            targets += ('</div><p>click on a cell to see the constant and the '
                'command-line invocation used to render it!<p>to recreate the '
                'entire image, use <p><code id="invocation">'
                + generateCLIinvocation(c, showcells=True)
                + '</code>')
        else:
            out.write('<img ')
            targets = (
                '<p><code id="invocation">'
                + generateCLIinvocation(c, showcells=True)
                + '</code>\n<p>click on the render to update the command with '
                'a new center (<code>-e</code>)'
            )
            # if there's only one cell, make clicking on the image update the
            # code#invocation tag with new coordinates centered on the click
            # location in the image
            end_script = ('\n<script>\n'
                f'let xmin = {graph["x"]["min"]},\n'
                f'dx = {graph["x"]["max"] - graph["x"]["min"]},\n'
                f'ymax = {graph["y"]["max"]},\n'
                f'dy = {graph["y"]["max"] - graph["y"]["min"]};')
            end_script += '''
(() => {
    let $ = e => document.getElementById(e),
        j = $('julia'),
        o = $('invocation');
    j.onclick = e => {
        let x = xmin + dx * e.offsetX / j.clientWidth;
        let y = ymax - dy * e.offsetY / j.clientHeight;
        o.innerHTML = o.innerHTML.replace(
            /-e (-?(\d|\.)+ ){2}/g, '-e '
            + x.toPrecision(4) + ' ' + y.toPrecision(4) + ' ');
    };
})();
</script>\n'''

        def tr(one, two):
            return f'<tr><td>{one}</td><td>{two}</td></tr>'

        # table with general render info
        out.write(
            f'src="../{fname}.png" id="julia">\n'
            + '<div id="container"><div id="contents">\n'
            + targets
            + '<table class="render-info">\n'
            + tr('z<sub>n + 1</sub>(z, c) =', orig_fn)
            + tr('c =', strcomplex(c))
            + tr('rendered area',
                f'({graph["x"]["min"]}, {graph["y"]["min"]}i) '
                f'to ({graph["x"]["max"]}, {graph["y"]["max"]}i)')
            + tr('center', f'({graph["x"]["c"]}, {graph["y"]["c"]}i)')
            + tr('zoom', f'{zoom}×')
            + tr('gradient speed', colorscale)
            + tr('escape radius', cutoff)
            + tr('iterations', iterations)
            + tr('c range', f'{crange:g}')
            + '</table>'
            + template_end.read()
            + end_script)

# if we're just regenerating the info, we're done
if no_render:
    exit()

# abstraction for ppm writing
def write_pixel(r, g, b, f):
    f.write(bytes([r, g, b]))

# formatting time-deltas in reasonable strings (why isn't this easier /
# built-in????)
def strtimedelta(time):
    return (f'{math.floor(time.seconds / 3600):02d}:'
        f'{math.floor((time.seconds % 3600) / 60):02d}:'
        f'{time.seconds % 60:02d}.{math.floor(time.microseconds / 1000):03d}')

with open(out_dir + '/' + fname + '.ppm', 'wb') as out:
    # magic numbers, etc
    out.write(bytes(f'P6\n{width} {height}\n255\n', encoding='ascii'))
    # initialize relevant variables
    z:      complex
    z_p:    complex
    color:  float
    graphy: float
    i:      int
    row:    int = 0
    col:    int = 0
    # for eta and total elapsed time estimation
    start = datetime.now()

    # iterate top to bottom
    for y in range(0, height):

        # if eta / % info requested, show it!
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

        # y component is the same for every x-coord, so we pre-calculate it
        graphy = stgY(y) * 1j
        # increment row if necessary, excepting moving past the edge of the
        # image
        if row != cellcount - 1 and y == yticks[row]:
            row += 1

        # we only *increment* the col index, so we have to reset it each row
        col = 0
        for x in range(0, width):
            if col != cellcount - 1 and x == xticks[col]:
                col += 1

            # z_0 is based on x and y coords
            # remember, stgX and stgY take columns and rows into account
            z_p = z = stgX(x) + graphy
            i = 0
            # smooth coloring using exponents????
            color = math.exp(-abs(z))
            # iterate, breaking if we exceed the orbit threshold
            for i in range(0, iterations):
                z = eval_fn(z, cgrid[col][row])
                if cmath.isnan(z) or not cmath.isfinite(z):
                    # oh no
                    # i can blame this on the user right?
                    z = z_p
                    break
                elif abs(z) > cutoff:
                    break
                z_p = z
                color += math.exp(-abs(z))

            # color is in the range of [0, iterations], so we have to
            # normalize it
            color /= iterations
            # we have way more magnitude info than can be expressed with 255
            # discrete values, so we oscillate light and dark as magnitude
            # increases, for a variety of colors. by scaling the blue
            # oscillator slightly longer and the red slightly shorter, we
            # get a nice rainbow!
            write_pixel(
                int(255 * math.sin(color * colorscale * 9 ) ** 2),
                int(255 * math.sin(color * colorscale * 10) ** 2),
                int(255 * math.sin(color * colorscale * 11) ** 2),
                out
            )

    # time elapsed
    end = datetime.now()
    print(f'Done! Completed in {strtimedelta(end - start)}')

if convert:
    # convert ppm to png (if requested, by default on)
    # don't load these libraries unless needed
    from subprocess import run
    print(f'Converting {fname}.ppm to {fname}.png')
    run(['magick', 'mogrify', '-format', 'png', f'{out_dir}/{fname}.ppm'])
    os.remove(f'{out_dir}/{fname}.ppm')

if open_html:
    # this works very smoothly, which is nice! thanks python!
    import webbrowser
    from urllib.request import pathname2url
    # only used once but boy oh boy does it make the code nicer
    def abspath(filename):
        return 'file:' + pathname2url(os.path.abspath(filename))

    webbrowser.open(abspath(f'{out_dir}/{info_dir}/{fname}.html'))
