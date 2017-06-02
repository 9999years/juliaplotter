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

## General Usage

```
usage: julia.py [-h] [--fn zₙ₊₁] [-c constant] [-a aspect] [-w width]
                [-i iterations] [-r c-range] [-e center center]
                [-n cell count] [-z zoom] [-g gradient speed] [-u escape]
```

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

### `--fn zₙ₊₁, -f zₙ₊₁`

The Julia set's function for iteration.

### `-c constant`

The constant c for the function zₙ₊₁(z, c). Enter `random` to select a random
value for c.

### `-a aspect, --aspect aspect`

The output image's w/h aspect ratio. Ex.: `-a 2` implies an image twice as wide
as it is tall.

### `-w width, --width width`

The output image's width.

### `-i iterations, --iterations iterations`

The iterations to calculate the set to.

### `-r c-range, --c-range c-range`

The range of c values to use — only relevant if the cell count option is used to
render a grid of sets; the c values for each sets will range from `(c_r -
crange, c_i - crange·i)` to `(c_r + crange, c_i + crange·i)`, where `c_r` and
`c_i` are the real and imaginary components of the constant supplied with `-c`.

### `-e cx cy, --center cx cy`

The coordinate the graph is centered around, entered as two floats separated by
a space. (Not a comma! No parenthesis!)

### `-n cell count, --cell-count cell count`

The number of rows and columns to render. A cell count of 1 will render a single
set, and other values will render grids of Julia sets. The different values of c
are determined by `--c-range` or `-r`.

### `-z zoom, --zoom zoom `

How zoomed in the render is. The distance between the center-point and the top /
bottom of the rendered area is 1 / zoom. Larger values of will produce a more
zoomed-in image, smaller values (<1) will produce a more zoomed-out image.

[1]: https://www.imagemagick.org/script/index.php
