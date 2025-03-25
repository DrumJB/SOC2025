import numpy as np

def create_inp_file(r, B, Te, pot, run_name):

	eps0 = 8.854e-12
	kB = 1.3806e-23
	e = 1.6021e-19
	n0 = 1e18
	amu = 1.6e-27
	mi = 2*amu
	qi = e
	cs = np.sqrt(kB*Te/mi)
	lamD = np.sqrt((eps0*kB*Te)/(n0*e**2))
	rL = (mi*cs)/(qi*B)
	if r/lamD < 1:
		print('!! Debyova delka vetsi nez polomer sondy !!')
	inp_file=f"""$plasma
	ksi={rL/lamD},
	tau=1,
	mu=200,
	P0=0.0
	PL=0.0
	alpha_xz=90.0,
	alpha_yz=90.0
$end

$geom
	Ly={int(10*rL/lamD)},
	Lz={int(10*rL/lamD)},
	dy=1
	dz=1
	Lz_low=0.0
	Lz_high=1000.0
	ta=3.10,
	tc=3.00,
	tp=4.10,
	Npc=50
	Npts_ratio=2.0
	scenario=2
	automatic_ta=0
$end

$control
	time_diag=.true.
	relax=.false.
	psolver=2
	ps_commsize=4
	diag_ntimes=10000
	dump_period=1000000
	pot_smoothing=0
	history_ntimes=10000
	ang_diag=.false.
	ang_bin=100.0
	potdebug=0
	pert_pot_file=''
$end

$optional
	param1=0.0
	param2=0.0
	param3=2.0
	param4=100.0
	param5=''
$end

$mks                                                                                                                                         
	mks_n0={n0}
	mks_Te={Te}
	mks_B={B}
	mks_main_ion_q=1.0
	mks_main_ion_m=2.0
	mks_par1=0.0
	mks_par2=0.0
	mks_par3=0.0                           
$end  

$domaindecomp
	no_slices_z=1,
	no_slices_y=1,
$end    

$num_blocks
	rectangles=0
	triangles=0
	circles=1
	shape1=0
	shape2=0
	shape3=0
$end

$circle
	name='pin'
	group=1
	ycentre={int(5*rL/lamD)},
	zcentre={int(5*rL/lamD)},
	radius={r/lamD},
	pot={pot}
	enable_erosion=0
	negative=.false.
	param1=0	
	param2=0
	param3=0
$end

$num_spec
	no_species=4
$end

$specie
	name='ions-top'
	T=1,
	m=1.0
	q=1.0
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'DF/fu_iz_TAU10.dat'
	injection_method = 0
	injection_rate_relative = .true.
	injection_rate = 0.5
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end
 
$specie
	name='electrons-top'
	T=1.0
	m=0.005
	q=-1.0
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'fu_ez.dat'
	injection_method = 0
	injection_rate_relative = .true.
	injection_rate = 0.5
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end 


$specie
	name='ions-bot'
	T=1,
	m=1.0
	q=1.0
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'DF/fu_iz_TAU10.dat'
	injection_method = 1
	injection_rate_relative = .true.
	injection_rate = 0.5
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end
 
$specie
	name='electrons-bot'
	T=1.0
	m=0.005
	q=-1.0
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'fu_ez.dat'
	injection_method = 1
	injection_rate_relative = .true.
	injection_rate = 0.5
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end 

$num_diag_regions
	no_diag_reg = 2
$end

$diag_reg
	y_low=1.0
	y_high=1.0
	z_low=1
	z_high=1
	diag_name = 'qnpot'
	record_property  = 1
	record_only_steady_state= .true.
	nbins = 1
	record_interval= 1.0
	diag_specie = 1
$end

$diag_reg
	y_low=1,
	y_high=3,
	z_low=2,
	z_high=4,
	diag_name = 'pin'
	record_property  = 1
	record_only_steady_state= .true.
	nbins = 1
	record_interval= 1.0
	diag_specie = 1
$end"""

	with open(f"{run_name}/r{r}B{B}Te{Te}Pot{pot}.inp", "w+") as output:
		output.write(inp_file)
		print('Inputs written.')

def create_sh_file(r, B, Te, pot, run_name):

	sh_file=f"""#!/bin/bash
#PBS -l nodes=1:soroban-node-02:ppn=1
#PBS -q long
#PBS -N r{r}B{B}Te{Te}Pot{pot}

HOST=$(hostname -s)

module purge
module load intelmpi/18
module load petsc/intel_oneapi
BIN_DIR=/compass/Shared/Common/IT/SW/SPICE/sw/soroban/unstable/var_te
EXE_NAME=spice-2.15-te-release.bin

RUN_CMD="mpirun -np 1"

ulimit -s unlimited

INP_DIR=/compass/home/buben/SPICE2/{run_name}
SCRATCH_DIR=/scratch/buben/{run_name}


DATA_DIR=$SCRATCH_DIR/d
BACKUP_DIR=$DATA_DIR-$RANDOM
INP_NAME=r{r}B{B}Te{Te}Pot{pot}
LOG_NAME=$SCRATCH_DIR/$INP_NAME.log

cd $BIN_DIR
mkdir -p $SCRATCH_DIR
mkdir -p $DATA_DIR
mkdir -p $BACKUP_DIR

cp -f $INP_DIR/$INP_NAME.inp $SCRATCH_DIR/$INP_NAME.inp
chmod a+rw $SCRATCH_DIR/$INP_NAME.inp

ldd $EXE_NAME

cp -vf $DATA_DIR/$INP_NAME*.mat $BACKUP_DIR
cp -vf $SCRATCH_DIR/$INP_NAME.log $BACKUP_DIR

echo "job $INP_NAME is running">$SCRATCH_DIR/$INP_NAME.running

$RUN_CMD $BIN_DIR/$EXE_NAME -v -i $SCRATCH_DIR/$INP_NAME.inp -t $DATA_DIR/$INP_NAME-t -o $DATA_DIR/$INP_NAME-o > $LOG_NAME

rm -f $SCRATCH_DIR/$INP_NAME.running"""

	with open(f"{run_name}/r{r}B{B}Te{Te}Pot{pot}.sh", "w+") as output:
		output.write(sh_file)

def create_launcher(rs, Bs, Tes, pots, run_name):

	with open(f"Launch-{run_name}.sh", "w+") as output:
		output.write('CHECK_MARK="✅"\n')
		output.write('echo -e "\n\e[4mRequest process\e[0m"\n')
		for r in rs:
			for B in Bs:
				for Te in Tes:
					for pot in pots:
						output.write(f"qsub {run_name}/r{r}B{B}Te{Te}Pot{pot}.sh\n")
						output.write('echo -e "${CHECK_MARK} '+str(f"r{r}B{B}Te{Te}Pot{pot}.sh")+'"\n')
		output.write('echo -e "All calculations have been requested."\n')
