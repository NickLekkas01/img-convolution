# Image-Convolution

## ΜPI - MPI & OPENMP

Applying image filters using Parallel Programming Techniques.

### Getting Started

Check your compiler's version of gcc by typing:

```gcc --version ```  - 4.4 or newer version is required.


Τo install mpicc. Type: 

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
``` mpiexec -f <machines> -n <processes> ./a.out width height type loops ```

##### MPI & OPENMP: 
``` mpiexec -n <processes> ./a.out width height type loops ``` 

### Εxamples

   Before                  |          After 50 loops
:-------------------------:|:-------------------------:
![](https://github.com/msiampou/Image-Convolution/blob/master/img/waterfall_grey.png)  |  ![](https://github.com/msiampou/Image-Convolution/blob/master/img/blur_wgrey.png)
![](https://github.com/msiampou/Image-Convolution/blob/master/img/waterfall.png)  |  ![](https://github.com/msiampou/Image-Convolution/blob/master/img/blur_w.png)

