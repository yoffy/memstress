# memory stress tool

## About This

This is a memory stressing tool. Used with commands such as `perf` and `time`.

## Requirements

* C++11 compiler

## Build

```
$ make
```

or

```
$ CXXFLAGS="-march=native" make
```

## Usage

```
$ ./memstress
```

## Examples

Compare cache misses for 524288 x 16KB (8GB) and 2048 x 4MB (8GB) workloads.
```
$ perf stat -ddd ./memstress $((2**19)) $(((2**4)*1024)) reverse8
$ perf stat -ddd ./memstress $((2**11)) $(((2**12)*1024)) reverse8
```
