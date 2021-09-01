# metagraphInterface
Interface class for Metagraph

### Usage

You need a build of `metagraph` (https://github.com/ratschlab/metagraph) somewhere on your system.

Copy the files in `include/mabl3` to your project. Use CMake for building to
save you from quite some pain. In your projects `CMakeLists.txt`, add

`add_subdirectory(include/mabl3)` (or wherever you have put the files)

This will take care of everything and make the `metagraphInterface` target
available for linking with your own targets.

When building your project, you will have to point CMake to the metagraph build,
see below for the option `-DMETAGRAPH_PROJECT_ROOT`.

## Build Graph

Use the script `script/createGraph.py` to create an annotated graph from Fasta files!

The interface requires that the metagraph was built and annotated in a certain way, above script will take care of doing this.
The script will create a `*.dbg` and a `*.column_coord.annodbg` file, those are the graph and annotation files you can use
with the `MetagraphInterface` class.

**Note:** The annotation only stores fasta sequence headers, thus they must be unique! `createGraph.py` will also output a
mapping of fasta file names and sequence headers in JSON format.

## Unit Tests

Clone this repo recursively

`$ git clone --recurse-submodules https://github.com/mabl3/metagraphInterface`

If you already cloned this repo but forgot the submodules, run

`$ git submodule update --init --recursive`

### Build

Out of source build using CMake

```
$ mkdir -p test/build && cd test/build
$ cmake .. -DMETAGRAPH_PROJECT_ROOT=/path/to/metagraph
$ make
```

Note that you need to build metagraph from source, this interface is
not compatible with the conda or docker images.

### Running Unit tests

#### Creating Testdata

CMake will do this step for you when compiling the unit tests

You can also do this manually (not necessary anymore when using CMake):

```
$ cd test/testdata
$ python3 generateTestdata.py --metagraph /path/to/metagraph/build/metagraph
```

This will create a bunch of files with sequences in them as well as an
annotated graph from these sequences.

#### Running Unit Tests

Then, in the build directory, run

`$ ./testMetagraphInterface`
