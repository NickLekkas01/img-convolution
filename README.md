# Image-Convolution

## ΜPI - MPI & OPENMP

Applying image filters using Parallel Programming Techniques.

### Getting Started

Check your compiler's version of gcc by typing:

```gcc --version ```  - 4.4 or newer version is required.


You also need to install. Type: 

```apt install mpich```

### Downloading

Download source code by typing:

``` git clone https://github.com/msiampou/Image-Convolution.git ```

### Compilation

##### MPI: 
``` mpicc mpi.c -lm ```

##### MPI & OPENMP: 
``` mpicc -fopenmp mpi.c -lm ```

### Running

##### MPI: 
``` mpiexec -f <machines> -n <processes> ./a.out width hight type loops ```

##### MPI & OPENMP: 
``` mpiexec -n <processes> ./a.out width hight type loops ``` 

