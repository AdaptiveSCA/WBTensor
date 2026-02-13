# WBTensor

A Tensor-Based Implementation for SIMON32/64.

Compiled with the support of [M4RI](https://bitbucket.org/malb/m4ri) for matrix computation.

# Install (Linux (Ubuntu/Debian))

```
$ sudo apt-get install libm4ri-dev

```

# Build

```
$ gcc main.c wb_simon.c -o wb_simon -lm4ri -std=c99
```

## Run

```
$ ./wbsimon
```