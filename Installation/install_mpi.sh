# 1) Download MPICH and unzip
wget https://www.mpich.org/static/downloads/4.0.2/mpich-4.0.2.tar.gz
tar xfz mpich-4.0.2.tar.gz

# 2) Make installation directory and temporary build directory
mkdir mpich-install
mkdir mpich-4.0.2-temp-build

# 3) Configure MPICH
export MAINDIR=$(pwd)
cd mpich-4.0.2-temp-build
../mpich-4.0.2/configure -prefix=$MAINDIR/mpich-install 2>&1 | tee c.txt

# 4) Build MPICH (takes at least 50 minutes)
time make 2>&1 | tee m.txt

# 5) Install MPICH commands
make install 2>&1 | tee mi.txt

# 6) Add to Install Path
export PATH=$MAINDIR/mpich-install/bin:$PATH

# 7) Create new make/local file to tell Torsten which MPI to use
cd ../Torsten/cmdstan
echo "TORSTEN_MPI=1" > make/local
echo "TBB_CXX_TYPE=gcc" >> make/local
echo "CXX=mpicxx" >> make/local
echo "CC=mpicc" >> make/local
echo "CXXFLAGS += -isystem $MAINDIR/mpich-install/include" >> make/local

# 8) Make Stan model with MPI
make clean-all
make ../example-models/twocpt_population/twocpt_population

# 9) Run
cd ../example-models/twocpt_population
mpiexec -n 2 ./twocpt_population \
  sample num_samples=50 num_warmup=50 algorithm=hmc engine=nuts max_depth=1 \
  data file=twocpt_population.data.R init=twocpt_population.init.R
  

# 10) Create hostfile and test on worker nodes
qconf -sel > hostfile

mpiexec -n 4 -bind-to core -f hostfile -l ./twocpt_population \
  sample num_samples=50 num_warmup=50 algorithm=hmc engine=nuts max_depth=1 \
  data file=twocpt_population.data.R init=twocpt_population.init.R
  