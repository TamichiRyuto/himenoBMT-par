# himenoBMT-par
This is a repository for parallel implimentation of himeno benchmark.

# Structure
```
.
├── Convergence
│   ├── himenoBMTconv-omp.cpp
│   ├── himenoBMTconv-sender.cpp
│   └── himenoBMTconv-stdpar.cpp
├── Normal
│   ├── himenoBMTxpa-omp.cpp
│   ├── himenoBMTxpa-sender.cpp
│   ├── himenoBMTxpa-stdpar
│   ├── himenoBMTxpa-stdpar.cpp
│   └── himenoBMTxpa.c
├── Makefile
└── README.md
```

`Normal/` directory contains C++ source files of himenoBMT based on the original implimentation.
`Normal/himenoBMTxpa.c` is the original version himenoBMT source code dynamically allocates arrays as pointer-based structures in the heap.
However, this source code shows performance in terms of MFLOPS, so we have corrected it to GFLOPS values.
`Convergence/` directory contains C++ source files of himenoBMT with convergence check.

### The file naming conventions
- `-omp` : The parallel implimentation using OpenMP
- `-stdpar` : The parallel implimentation using C++ Standard Parallelism
- `-sender` : The asynchronus parallel implimentation using C++ Sender/Receiver

Sender/Receiver implimentations support only NVIDIA GPUs. 
See [NVIDIA/stdexec GitHub page](https://github.com/NVIDIA/stdexec/) for more information.

# How to Benchmark
Assume a docker-enabled Linux environment.
```
$ git clone https://github.com/TamichiRyuto/himenoBMT-par
```

### Launch the container with the `docker run` command
```
$ docker run --gpus all -it --name himenoBMT-par -v $(pwd)/himenoBMT-par:/himenoBMT-par --workdir /himenoBMT-par nvcr.io/nvidia/nvhpc:24.9-devel-cuda12.6-ubuntu24.04
```

### Compile source files
```
$ make -j4
```

### Running benchmark programs
Problem size is selectable follow.

| Problem size | Array size |
|:-------------:|:----------|
| XS | $(64 \times 32 \times 32)$ |
| S | $(128 \times 64 \times 64)$ |
| M | $(256 \times 128 \times 128)$ |
| L | $(512 \times 256 \times 256)$ |
| XL | $(1024 \times 512 \times 512)$ |

#### Original himenoBMT
```
$ ./himenoBMTxpa <size>
$ ./himenoBMTxpa-omp <size>
$ ./himenoBMTxpa-stdpar <size>
$ ./himenoBMTxpa-sender <size>
```

#### himenoBMT with convergence check
```
$ ./himenoBMTconv-omp <size>
$ ./himenoBMTconv-stdpar <size>
$ ./himenoBMTconv-sender <size>
```
