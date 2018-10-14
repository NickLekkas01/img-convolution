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

### Εxamples

<p float="left">
  <img src="https://github.com/msiampou/Image-Convolution/blob/master/img/waterfall_grey_1920_2520.raw" width="100" />
  <img src="https://github.com/msiampou/Image-Convolution/blob/master/img/blur_wgrey.png" width="100" /> 
</p>

   Before                  |          After
:-------------------------:|:-------------------------:
![](https://github.com/msiampou/Image-Convolution/blob/master/img/waterfall_grey_1920_2520.raw)  |  ![](https://github.com/msiampou/Image-Convolution/blob/master/img/blur_wgrey.png)
