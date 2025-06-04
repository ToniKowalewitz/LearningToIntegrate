# Unit tests and regression tests #


## Running tests locally ##

Type

```matlab
runtests('tests', 'IncludeSubPackages', true)
```

or selectively

```matlab
runtests('tests.unit')
runtests('tests.demo')
runtests('tests.style')
runtests('tests.unit.TestAssembleLaplaceSensitivity')
```

from matlab interpreter, or

```shell
matlab -nodisplay -r run_tests
```

from shell. Always work in repository root dir.


## CI ##

[![pipeline status](https://git.tu-freiberg.de/ng/curlcurl/badges/master/pipeline.svg)](https://git.tu-freiberg.de/ng/curlcurl/commits/master)

Tests are also run automatically on every push to reposity by
CI service <https://git.tu-freiberg.de/ng/curlcurl/pipelines>.
Commit sitting on tip of the branch is tested. New commits
being parents of the tip are not.

To prevent running CI (because it is known that the commit
is broken or it is just a change in docs/comments which does
not need testing) drop the phrase `[skip ci]` in the commit
message on the tip of the branch to be pushed.
