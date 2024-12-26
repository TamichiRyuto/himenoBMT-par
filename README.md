# himenoBMT-par
This is a repository for parallel implimentation of himeno benchmark.

# Structure
**.**   
**├── Convergence** *#himenoBMT C++ sources with convergence check*  
**│   ├── himenoBMTconv-omp.cpp** *#OpenMP implimentation*   
**│   ├── himenoBMTconv-sender.cpp** *#Sender/Receiver implimentation (NVIDIA's GPU Only)*  
**│   └── himenoBMTconv-stdpar.cpp** *#Stdpar implimentation*  
**├── Normal** *#himenoBMT C++ sources based on the original implimentation*  
**│   ├── himenoBMTxpa-omp.cpp** *#OpenMP implimentation*  
**│   ├── himenoBMTxpa-sender.cpp** *##Sender/Receiver implimentation (NVIDIA's GPU Only)*  
**│   ├── himenoBMTxpa-stdpar.cpp** *#Stdpar implimentation*  
**│   └── himenoBMTxpa.c** *#Original himenoBMT C language source*  
**└── README.md**  
