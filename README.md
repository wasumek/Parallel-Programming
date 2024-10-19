# Parallel Programming Projects

This repository contains projects developed as part of a course on parallel programming. The projects focus on leveraging multiple processors or cores to solve computationally intensive problems more efficiently. This is achieved using different parallel programming models and frameworks such as OpenMP and MPI (Message Passing Interface).

## Table of Contents
- [Overview](#overview)
- [Projects](#projects)
- [Installation](#installation)
- [Usage](#usage)

## Overview

Parallel programming is a paradigm that allows the simultaneous execution of multiple processes or threads to improve computational performance, especially in tasks that are computationally expensive. The use of parallel architectures has become fundamental in modern computing, from high-performance scientific applications to real-time systems. This repository includes implementations of various parallel algorithms and explores different parallel programming models, with a focus on improving computation speed and efficiency on multicore systems.

## Projects

1. **OpenMP for Multithreading**:
    - This project demonstrates how OpenMP can be used to parallelize code at the thread level. By using directives in C or C++, loops and other code sections are parallelized, showcasing the performance improvements achievable with multicore CPUs. The project includes examples of matrix multiplication and numerical simulations.

2. **Message Passing Interface (MPI) for Distributed Memory Systems**:
    - MPI is utilized in this project to parallelize tasks across multiple processors in a distributed memory environment. MPI allows processes to communicate by passing messages, making it suitable for large-scale computations. Examples include solving large linear systems and parallel matrix operations.

## Installation

To run the code in this repository, clone it to your local machine:
```bash
git clone https://github.com/wasumek/Parallel-Programming.git
```

### Dependencies:
- OpenMP (for multithreading projects)
- MPI (for distributed memory parallelism)

## Usage

Follow the instructions in the respective project folders to compile and execute the parallelized code. For OpenMP and MPI projects, you can compile the C/C++ code with GCC or other compilers that support parallel frameworks:
```bash
gcc -fopenmp -o outputfile sourcefile.c   # For OpenMP
mpicc -o outputfile sourcefile.c          # For MPI
```
