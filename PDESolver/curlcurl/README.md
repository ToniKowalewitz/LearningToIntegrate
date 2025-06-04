# (not only) curl-curl solver #


## Installation ##

In matlab interpreter, just `cd` to repository root directory.
You need Matlab with Symbolic Math Toolbox. Supported versions
are listed in [CI configuration file](.gitlab-ci.yml).


## Quick start ##

Crank up matlab, `cd` to repository root and type `demo.cavity`.


## Documentation ##

Documentation is built on CI sevice on every push, deployed and
hosted on <http://ng.git-pages.tu-freiberg.de/curlcurl/>.

You can build the documentation locally by running
`matlab -nodisplay -r build_doc` from shell. For that you need
to install [patched version of m2html](https://github.com/blechta/m2html)
by, e.g.,

```shell
git clone https://github.com/blechta/m2html.git
matlab -r "addpath(strcat(pwd, '/m2html')); savepath;"
```

## Authors ##

* Jan Blechta <jan.blechta@math.tu-chemnitz.de>
* Felix Eckhofer <felix@eckhofer.com>
* Toni Kowalewitz <toni.kowalewitz@mathematik.tu-chemnitz.de>
* Mathias Scheunert <mathias.scheunert@math.tu-freiberg.de>


## License ##

MIT License
```
Copyright (C) 2018-2020 curl-curl authors

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```


Note that the curl-curl source code contains code from other projects
in `+util` licensed as BSD 2-clause. Such files contains explicit
copyright/license header.


## CI ##

* tests [![pipeline status](https://gitlab.hrz.tu-chemnitz.de/geosax/curlcurl/badges/master/pipeline.svg)](https://gitlab.hrz.tu-chemnitz.de/geosax/curlcurl/commits/master)
* docs deploy [![pipeline status](https://git.tu-freiberg.de/ng/curlcurl/badges/master/pipeline.svg)](https://git.tu-freiberg.de/ng/curlcurl/commits/master)
