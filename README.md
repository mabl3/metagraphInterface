# metagraphInterface
Interface class for Metagraph

### Build

Out of source build using CMake

```
$ mkdir -p test/build && cd test/build
$ cmake .. -DMETAGRAPH_PROJECT_ROOT=/path/to/projects2014-metagenome
$ make
```

If you built `metagraph` with Folly but Folly is in a non-standard location,
specify the path to `libfolly.a` in the `cmake` command using `-DFOLLY_PATH=/path/to/folly/lib`

### Running Unit tests

In build directory

`$ ./testMetagraphInterface`
