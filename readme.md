# Python Julia Plotter

A script to render Julia sets or grids of Julia sets with a comfortable
command-line interface, enabled-by-default conversion to `.png` ([`magick`][1]
required in path), and pretty-by-default HTML export.

## Why?

* **Low memory requirements** — Instead of storing the image data  an array
  before writing it, it’s written one byte at a time. This means that rendering
  a 65536×65536 image takes up just as much (run-time) memory as rendering a
  500×500 image.
* **Flexibility** — Any formula is accepted, with fairly versatile equation
  parsing.

```
usage: julia.py [-h] [--fn zₙ₊₁] [-c constant] [-a aspect] [-w width]
                [-i iterations] [-r c-range] [-e center center]
                [-n cell count] [-z zoom] [-g gradient speed] [-u escape]
```

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
