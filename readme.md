# Python Julia Plotter

A script to render [Julia sets](https://en.m.wikipedia.org/wiki/Julia_set) (the
sets of points in the complex plane that don’t escape with the repeated
application of a complex rational function) or grids of Julia sets with a
comfortable command-line interface, enabled-by-default conversion to `.png`
([`magick`][1] required in path), and pretty-by-default HTML export.

## Why?

### Low memory requirements

Instead of storing the image data  an array before writing it, it’s written one
byte at a time. This means that rendering a 65536×65536 image takes up just as
much (run-time) memory as rendering a 500×500 image. The raw `.ppm` output will
take up `9 + ⌈log10(w)⌉ + ⌈log10(⌊w/a⌋)⌉ + 3w⌊w/a⌋` bytes (three per pixel plus
a header, where `w` is the image width and `a` is the aspect ratio).


### Flexibility

Any formula is accepted, with fairly versatile equation parsing.

Accepted features:

* Imaginary numbers in the form of (ex. `1i`, `3.2i`, `.3i`, etc.)
* Coefficients for variables and functions (`2tan(z)`, `0.3c`, `12 pi`)
* Implicit parenthesis (ONLY if the variable passed to a function is unmodified;
  `tanz` and `log10 z` work, `sin z^2` will parse as `(sin z)^2`)
* Implicit multiplication (`(z - 2)(z + c)`, `z(z - 3)`, `2c z`, `sinz tanz`)
* Exponentiation (`z^2`, `z^-1`)

### Grid rendering

`julia.py` contains a built-in system for rendering *grids* of Julia sets. Why?
Because the same equation that makes this very interesting image:

![Julia set for f(z) = (z-c)(z+2-0.5i)(z+0.5c)(z), c =
0.9i](https://i.imgur.com/4QQRShw.png)

    ./julia.py -f "(z-c)(z+2-0.5i)(z+0.5c)(z)" -c "0 + 0.9 i" -i 32 -w 2048 -a 1.0 -e 0.0 0.0 -z 1.0 -g 1.0 -u 30.0

(Note: Most of those arguments are superfluous, but are outputted in case
defaults change. The actual rendering command was something much closer to
`./julia.py -f "(z-c)(z+2-0.5i)(z+0.5c)(z)" -c "0 + 0.9 i" -w 2048`)

Creates this very uninteresting image:

![Julia set for f(z) = (z-c)(z+2-0.5i)(z+0.5c)(z), c = 0.3 +
0.3i](https://i.imgur.com/76RVZby.png)

    ./julia.py -f "(z-c)(z+2-0.5i)(z+0.5c)(z)" -c "0.3 + 0.3 i" -i 32 -w 2048 -a 1.0 -e 0.0 0.0 -z 0.75 -g 1.0 -u 30.0

These two sets are in the same *family* of equations, but have a different
constant *c*. With the `-n cells` option, `julia.py` will render a grid of
cells × cells sets with c-values ranging from -range to range in the real axis,
and -range·i to range·i in the imaginary axis, where the range is controlled
with the `-r range` option. We may discover which constants are interesting and
which are not with a preliminary render:

![Grid of Julia sets for f(z) = (z-c)(z+2-0.5i)(z+0.5c)(z)](https://i.imgur.com/fs6Fuv6.png)

    julia.py -f "(z-c)(z+2-0.5i)(z+0.5c)(z)" -n 11 -i 32 -w 2048 -a 1.0
    -e 0 0 -z 0.5 -g 1.0 -u 30

(Aside: Although it might look it at first, the grid is not symmetrical over the
real or imaginary axis.)

`julia.py` will then automatically open the HTML output, allowing us to click on
an interesting cell and copy the command-line invocation to create a more
detailed render.

![Screenshot of HTML output](https://i.imgur.com/DmS6fex.png)

## HTML Output Features

At a glance, the HTML output contains everything you need to recreate a render,
the easy way (by copying and pasting the given command-line invocation) or the
hard way (by manually typing in the arguments in the table [do not do this]).
There are two other notable features:

1. In renders with multiple cells (i.e. with `-n` > 1), clicking on a cell in
   the image will reveal the column, row, c-value, and a command-line invocation
   to render that cell in a larger image (by default to the same width as the
   render it came from). Note that this uses an image `<map>`, and is thus
   pixel-dependent; if the render is larger than your screen, it will be resized
   down, and Javascript will have to be enabled for the image map to resize
   correctly (easy/possible thanks to David J. Bradshaw’s
   [`imageMapResize.js`](https://github.com/davidjbradshaw/image-map-resizer/))
2. In renders with one cell (i.e. with `-n 1` or with `-n` omitted), clicking on
   the render will update the shown command-line invocation with an updated
   center value (`-e`). By clicking on a region of the render you would like to
   see enlarged and increasing the zoom value (`-z`), more detailed renders can
   easily be created.

## General Usage

```
usage: julia.py [-h] [-f zₙ₊₁] [-c constant] [-a aspect] [-w width]
                [-i iterations] [-r c-range] [-e center center]
                [-n cell count] [-z zoom] [-g gradient speed] [-u escape]
```

Note that if `-f` is omitted, it will default to `z^2 + c`, and if `-c` is
omitted, it will default to `0 + 0i`, but if *both* `-f` and `-c` are omitted,
`-f` will default to `z^2 + c` and **`-c` will default to `random`** — this may
be unexpected, but results in interesting-and-unpredictable-by-default renders.

## Example Renders

![Julia set for f(z) = (z^3)/(((z-2)^5)(z+0.5i)) + c, c = 7/6 +
1/6i](https://i.imgur.com/dEbYTN8.png)

    ./julia.py -f "(z^3)/(((z-2)^5)(z+0.5i)) + c" -c " 1.16667 + 0.166667i" -i 32 -w 2048 -a 1.0 -e 2.0 0.0 -z 0.5 -g 1.0 -u 30

![Julia set for f(z) = z^2 + c, c = 0.285 + 0.01i](https://i.imgur.com/So2smbd.png)

    ./julia.py -c "0.285 + 0.01i" -z 0.75 -w 2048

## Arguments and options

The following is an enumeration of most of the useful and common arguments.
There are a few others and this readme may fall slightly out of date — view the
most up-to-date and complete list with `-h` or `--help`.

### Quick argument table

Arg   | Long arg       | What it controls                | Default value
----- | ----------     | ------------------              | --------------
`-f`  | `--fn`         | Render function zₙ₊₁(z, c)      | `z^2 + c`
`-c`  | `--constant`   | Render constant (`c` in `-f`)   | `0 + 0i`
`-a`  | `--aspect`     | Image aspect ratio              | `1.0`
`-w`  | `--width`      | Image width                     | `500`
`-i`  | `--iterations` | Fractal iterations              | `32`
`-r`  | `--c-range`    | Range of c-values when `-n` > 1 | `1.5`
`-n`  | `--cell-count` | Cell count                      | `1`
`-e`  | `--center`     | Render center                   | `0 0`
`-z`  | `--zoom`       | Image zoom                      | `1.0`
`-s`  | `--silent`     | Info output, shelling out       | Off

### `--fn zₙ₊₁, -f zₙ₊₁` (Default: `z^2 + c`)

The Julia set's function for iteration.

### `-c constant` (Default: `0 + 0i`)

The constant c for the function zₙ₊₁(z, c). Enter `random` to select a random
value for c.

### `-a aspect, --aspect aspect` (Default: `1`)

The output image's w/h aspect ratio. Ex.: `-a 2` implies an image twice as wide
as it is tall.

### `-w width, --width width` (Default: `500`)

The output image's width.

### `-i iterations, --iterations iterations` (Default: `32`)

The iterations to calculate the set to.

### `-r c-range, --c-range c-range` (Default: `1.5`)

The range of c values to use — only relevant if the cell count option is used to
render a grid of sets; the c values for each sets will range from `(c_r -
crange, c_i - crange·i)` to `(c_r + crange, c_i + crange·i)`, where `c_r` and
`c_i` are the real and imaginary components of the constant supplied with `-c`.

### `-n cell count, --cell-count cell count` (Default: `1`)

The number of rows and columns to render. A cell count of 1 will render a single
set, and other values will render grids of Julia sets. The different values of c
are determined by `--c-range` or `-r`.

### `-e cx cy, --center cx cy` (Default: `0 0`)

The coordinate the graph is centered around, entered as two floats separated by
a space. (Not a comma! No parenthesis! Technically two separate arguments
consumed by one option.)

### `-z zoom, --zoom zoom` (Default: `1`)

How zoomed in the render is. The distance between the center-point and the top /
bottom of the rendered area is 1 / zoom. Larger values of will produce a more
zoomed-in image, smaller values (<1) will produce a more zoomed-out image.

### `-s, --silent` (Default: Off)

Don't log info, show progress, convert the .ppm to a .png, or open the file when
finished. Equivalent to `--no-open --no-convert --no-progress --no-info`.

## License

MIT / Expat; see [license.txt](blob/master/license.txt).

[1]: https://www.imagemagick.org/script/index.php
