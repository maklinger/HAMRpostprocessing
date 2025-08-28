# To get H-AMR running on HIPSTER

First, login to HIPSTER from within an UvA network or VPN (with your uvanetID):
```shell
ssh UvanetID@hipster.science.uva.nl
```

## 1. H-AMR makefile
The Makefile that worked for me with the RHAMR_CUDA_CIP_OCL branch is:
```make
#### set USEICC to 0 if you want gcc compiler options, else set to 1 to use icc
########  gcc generally used for debugging with -g option so we can use gdb 
USEICC = 0
ifeq ($(USEICC),0)
CC       = CC 
CCFLAGS  = -I${ROCM_PATH}/include -std=c++11 -D CRAY_CPU_TARGET=x86-64 -D__HIP_ROCclr__ -D__HIP_ARCH_GFX942__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx942 -fgpu-rdc -x hip -O3
endif

EXTRALIBS = -lm -fgpu-rdc --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64 -lstdc++
CC_COMPILE  = $(CC) $(CCFLAGS) -c 
CC_LOAD     = $(CC) $(CCFLAGS)


CUDA_COMPILE  = hipcc -fgpu-rdc --offload-arch=gfx942 -I${MPICH_DIR}/include -I${ROCM_PATH}/include -D__HIP_PLATFORM_AMD__ -fopenmp -O3 -c
CUDA_LOAD  = hipcc -v -fgpu-rdc --offload-arch=gfx942 -L${MPICH_DIR}/lib -lmpi -L${MPICH_DIR}/lib -lmpi -fopenmp
GPU_FILES = GPU_boundcomP.cu GPU_boundcomF.cu GPU_boundcomE.cu GPU_main.cu GPU_program1.cu GPU_program2.cu
.c.o:
	$(CUDA_COMPILE) $*.c

EXE = harm
all: $(EXE)
    
OBJS = \
AMR.o boundcomB.o boundcomE.o boundcomF.o boundcomP.o balance_load.o\
bounds.o coord.o const_trans.o const_trans_res.o diag.o dump.o eos_helm.o fixup.o\
GPU_boundcomE.o GPU_boundcomP.o GPU_boundcomF.o GPU_program1.o GPU_program2.o GPU_main.o\
hllc.o LAS.o init_collapsar.o init_mag.o init_misc.o init_nsm.o init_tde.o \
init_tests.o init_thindisk.o init_torus.o init_torus_grb.o init_torus_spherical.o \
interp.o lu.o main.o memory.o metric.o phys_2T.o phys_mhd.o phys_neutrinos.o \
phys_nu.o phys_nuclear.o phys_radiation.o phys_res.o radiation.o ranc.o restart.o step_ch.o step_ch_res.o \
u2p_util.o utoprim_1dfix1.o utoprim_1dvsq2fix1.o utoprim_2d.o utoprim_3d_res.o utoprim_nm.o wrapper.o

INCS = \
decs.h decs_MPI.h decsCUDA.h defs.h include.h u2p_defs.h u2p_util.h config.h
$(OBJS) : $(INCS) makefile
$(EXE): $(OBJS) $(INCS) makefile
	$(CUDA_COMPILE) $(GPU_FILES)
	$(CUDA_LOAD) $(OBJS) -o $(EXE)
clean:
	/bin/rm -f *.o *.il
	/bin/rm -f $(EXE) image_interp

```
Note that I removed the part to print the user summary.


## 2. Running H-AMR
This can be used in two ways: via a script or interactively on the fp64 node to run/debug simulations on the GPUs

The HIPSTER architecture makes running H-AMR a bit complex. The logic is, that you first submit a job to slurm to get access to the fp64 node, and then start the apptainer (I used the one from Jeroen Roodhart), in which you then compile and run the H-AMR code.
Jeroen's apptainer is located at:
`/opt/sw-fnwi/sif/pytorch_rocm6.2.3_ubuntu22.04_py3.10_pytorch_release_2.3.0.sif`


### A. Running H-AMR on fp64 via a script

For readability, I generated two job files in the same directory as the rest of the H-AMR code. The first one contains the slurm part and starts the second one inside the apptainer instance. First, in `run_job.sh`:

```shell
#!/bin/bash
#SBATCH --job-name=test
#SBATCH --nodes=1               # Number of nodes - leave at 1
#SBATCH --ntasks=6              # Number of tasks - leave at 6 (AMD MI300A has 6 GPU-like units)
#SBATCH --cpus-per-task=24      # Number of cpus - leave at 24 (AMD MI300A has 24 CPUs)
#SBATCH --output=/home/....../output.out   # Output file (output to home dir so you can see progress)
#SBATCH --error=/home/......./output.err    # Error file (output to home dir so you can see progress)
#SBATCH --partition=fp64        # Check partition with `sinfo` - leave at fp64
#SBATCH --gres=gpu:2            # Request 2 GPUs - adjust to at most 4 (leaves no access to GPU anymore to inspect first results while code is still running)
#SBATCH --time=1-00:00:00       # Max runtime: eg. 1 day
#SBATCH --mem=256G              # Memory limit (128 per GPU)


# Define paths (no "/" at end) - adjust to your case
UNAME="mklinge1"
HOME_PATH="/home/${UNAME}/HAMR"
SETUP_NAME="RHAMR_CUDA_CIP_OCL_blastwave"
SIF_PATH="/opt/sw-fnwi/sif/pytorch_rocm6.2.3_ubuntu22.04_py3.10_pytorch_release_2.3.0.sif"
SCRATCH_DIR="/local_scratch"

echo "ntasks=" ${SLURM_NTASKS:-1} "CPUs/task=" ${SLURM_CPUS_PER_TASK:-1}

# Run commands inside Apptainer
# Use srun from the host to launch into the container
apptainer exec \
  --bind $HOME \
  --bind $SCRATCH_DIR:$SCRATCH_DIR \
  --rocm $SIF_PATH \
  "${HOME_PATH}/${SETUP_NAME}/run_in_container.sh" ${SLURM_NTASKS:-1} ${SLURM_CPUS_PER_TASK:-1} \
  $UNAME $HOME_PATH $SETUP_NAME $SCRATCH_DIR

```
It is important to bind any path you want to access from within the apptainer instance and to set the `--rocm` flag to enable the ROCM compiler and GPU functionality.

The second script `run_in_container.sh` contains:
```shell
#!/bin/bash
export MPICH_DIR=/opt/ompi
echo 'Using MPICH_DIR=' $MPICH_DIR
export OMP_NUM_THREADS=$2
echo 'Using OMP_NUM_THREADS=' $2

# Copy files to local scratch folder
UNAME=$3
HOME_PATH=$4
SETUP_NAME=$5
SCRATCH_DIR=$6
cp -r "${HOME_PATH}/${SETUP_NAME}/" "${SCRATCH_DIR}/${UNAME}/${SETUP_NAME}"
echo "Copied!"

# Go to local scratch folder
cd "${SCRATCH_DIR}/${UNAME}/${SETUP_NAME}"

# Wait for a second to make sure everything is copied
sleep 1s

# Compile code
echo "Start compilation.."
make clean
make -j
echo 'compiled successfully!'

# This needs the summary script and is optional
# make summary 

# Run the code with mpi
echo 'Executing mpirun -n ' $1 ' ./harm '
mpirun -n $1 ./harm 

```

### B. Running H-AMR on fp64 interactively

This is a bit more straight forward:
```console
srun -n1 --gres=gpu:1 -p fp64 --mem=128G --time=01:30:00 --pty apptainer shell --bind $HOME --bin /local_scratch:/local_scratch --rocm /opt/sw-fnwi/sif/pytorch_rocm6.2.3_ubuntu22.04_py3.10_pytorch_release_2.3.0.sif
```

Then you can go to the H-AMR directory and compile and run:
```console
make clean
make -j
mpirun -n4 ./harm
```

## 3. Postprocessing

### Copy to mid-term storage
This way of running H-AMR just runs it until the job's time is over (or the memory is used up...), such that the result files are still on the fp64 nodes' local scratch. To get access to that node and copy them back, I use the follwoing alias (basically an interactive run):

```shell
interactivefp64(){
    srun -n1 -p fp64 --mem=10G --time=02:30:00 --gres=gpu:1 --pty apptainer shell \
        --bind $HOME --bind /local_scratch:/local_scratch \
        --bind /scratch:/scratch --bind /fnwi_fs/:/fnwi_fs/ \
        --rocm /opt/sw-fnwi/sif/pytorch_rocm6.2.3_ubuntu22.04_py3.10_pytorch_release_2.3.0.sif
}
```
Note that you have to request at least one GPU for a job, so this only works when the cluster is not busy. It might better at some point to trigger some automatic copying..

After you ran your simulation, you can copy it to: `/fnwi_fs/fnwi/hipster`. Just make your own folder and keep in mind that we only have about 10TB (per user or per group..).


### Create an apptainer for postprocessing
For simple plotting/postprocessing on the lgin node, I created my own apptainer. For this, I used the following file `ubuntu22_python.def` just somewhere in a separate folder:

```
Bootstrap: docker
From: ubuntu:22.04

%environment
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=/opt/conda/envs/analysis/bin:$PATH
    export CONDA_PREFIX=/opt/conda/envs/analysis
    export JUPYTERLAB_PORT=12345

%post
    apt-get update && apt-get install -y \
        curl \
        wget \
        git \
        build-essential \
        ca-certificates \
        gcc \
        g++ \
        make \
        bash \
        bzip2 \
        libssl-dev \
        htop \
        && apt-get clean

    # Download micromamba
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

    # Create conda environment at specific prefix (no -n!)
    micromamba create -y -p /opt/conda/envs/analysis python=3.11 \
        numpy \
        scipy \
        sympy \
        matplotlib \
        pandas \
        xarray \
        h5py \
        ipython \
        jupyterlab \
        ipympl \
        ipywidgets \
        pip \
        numba \
        cython \
        compilers \
        setuptools \
        wheel \
        && micromamba clean -a -y

%runscript
    #!/bin/bash
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=/opt/conda/envs/analysis/bin:$PATH
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate analysis
    exec "$@"

%startscript
    #!/bin/bash
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH=/opt/conda/envs/analysis/bin:$PATH
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate analysis
    exec jupyter-lab --no-browser --ip 0.0.0.0 --port $JUPYTERLAB_PORT
```

This creates a container with micromamba and a python environment with the packages listed. You can add packages and recompile. Compiling the container works as:
```console
apptainer build ubuntu22_python.sif ubuntu22_python.def
```

This can be used interactively in VS code by starting/stopping an apptainer instance to keep running.

Start the server like: 
jupyterserverstart() {
    > /home/mklinge1/.apptainer/instances/logs/$(hostname)/$USER/jupyterlab.err
    apptainer instance start --bind /fnwi_fs/:/fnwi_fs/ /home/mklinge1/ApptainerImages/ubuntu22_python.sif jupyterlab
    sleep 5s
    cat /home/mklinge1/.apptainer/instances/logs/$(hostname)/$USER/jupyterlab.err | grep -o 'http.*lab?token.*' | sort | uniq
}

Then copy the printed link and use that in VS code. When done, stop the instace via:

```shell 
jupyterserverstop() {
    apptainer instance stop jupyterlab
}
```