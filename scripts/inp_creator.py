#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 18:49:05 2023

@author: jachym
"""
import numpy
def write_new_file(type, tau, H, He, pot, name="gamma-TPC1223"):
    eps0 = 8.854e-12
    kB = 1.3806e-23
    e = 1.6021e-19
    Te = 20*e/kB
    n0 = 1e18
    B = 2.5
    amu = 1.6e-27
    mi = round(4.0*int(He)/100+1.0*int(H)/100,2)*amu
    qi = round(2.0*int(He)/100+1.0*int(H)/100,2)*e
    cs = numpy.sqrt(Te*kB/mi)
    lamD = numpy.sqrt((eps0*kB*Te)/(n0*e**2))
    rL = (mi*cs)/(qi*B)
    inp_file01=f"""$plasma
	ksi={round(rL/lamD,2)},
	tau={int(tau)/10},
	mu=200,
	P0=-{int(pot)}.0d0
	PL=0.0
	alpha_xz=90,
	alpha_yz=30.0,
$end

$geom
	Ly=30,
	Lz=80,
	dy=0.50
	dz=0.50
	Lz_low=0.0
	Lz_high=1000.0
	tc=3.00,
	ta=3.50,
	tp=4.00,
	Npc=50
	Npts_ratio=2.0
	scenario=1
	automatic_ta=0
$end

$control
	time_diag=.true.
	relax=.false.
	psolver=2
	ps_commsize=4
	diag_ntimes=50000
	dump_period=50000
	pot_smoothing=0
	history_ntimes=50000
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
	mks_n0=1.00E+18
	mks_Te=20
	mks_B=2.5
	mks_main_ion_q={round(2.0*int(He)/100+1.0*int(H)/100,2)}
	mks_main_ion_m={round(4.0*int(He)/100+1.0*int(H)/100,2)}
	mks_par1=0.0
	mks_par2=0.0
	mks_par3=0.0
$end

$domaindecomp
	no_slices_z=1,
	no_slices_y=2,
$end

$num_blocks
	rectangles=1
	triangles=0
	circles=0
	shape1=0
	shape2=0
	shape3=0
$end

$rectangle
	name='body'
	group=1
	ylow=0,
	zlow=0,
	yhigh=30,
	zhigh=15,
	pot=-{int(pot)}.0
	enable_erosion=0
	negative=.false.
	param1=0
	param2=0
	param3=0
$end

$num_spec
	no_species=2
$end

$specie
	name='ions-base-H{H}He{He}'
	T=1,
	m={round(2.0*int(He)/100+0.5*int(H)/100,2)}
	q={round(2.0*int(He)/100+1.0*int(H)/100,2)}
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'DF/fu_iz_TAU{tau}.dat'
	injection_method = 0
	injection_rate_relative = .true.
	injection_rate = 1.0
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end

$specie
	name='electrons-base'
	T=1.0
	m=0.005000000
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

$num_diag_regions
	no_diag_reg = 1
$end

$diag_reg
	y_low=14.0
	y_high=16.0
	z_low=13
	z_high=15
	diag_name = 'pot'
	record_property  = 1
	record_only_steady_state= .false.
	nbins = 1
	record_interval= 1.0
	diag_specie = 1
$end

"""

    mi = 4.0*amu
    qi = 2.0*e
    cs = numpy.sqrt(Te*kB/mi)
    lamD = numpy.sqrt((eps0*kB*Te)/(n0*e**2))
    rL = (mi*cs)/(qi*B)
    inp_file02=f"""$plasma
	ksi={round(rL/lamD,2)},
	tau={int(tau)/10},
	mu=200,
	P0=-{int(pot)}.0d0
	PL=0.0
	alpha_xz=90,
	alpha_yz=30.0,
$end

$geom
	Ly=30,
	Lz=80,
	dy=0.50
	dz=0.50
	Lz_low=0.0
	Lz_high=1000.0
	tc=3.00,
	ta=3.50,
	tp=4.00,
	Npc=50
	Npts_ratio=2.0
	scenario=1
	automatic_ta=0
$end

$control
	time_diag=.true.
	relax=.false.
	psolver=2
	ps_commsize=4
	diag_ntimes=50000
	dump_period=50000
	pot_smoothing=0
	history_ntimes=50000
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
	mks_n0=1.00E+18
	mks_Te=20
	mks_B=2.5
	mks_main_ion_q=2.0
	mks_main_ion_m=4.0
	mks_par1=0.0
	mks_par2=0.0
	mks_par3=0.0
$end

$domaindecomp
	no_slices_z=1,
	no_slices_y=2,
$end

$num_blocks
	rectangles=1
	triangles=0
	circles=0
	shape1=0
	shape2=0
	shape3=0
$end

$rectangle
	name='body'
	group=1
	ylow=0,
	zlow=0,
	yhigh=30,
	zhigh=15,
	pot=-{int(pot)}.0
	enable_erosion=0
	negative=.false.
	param1=0
	param2=0
	param3=0
$end

$num_spec
	no_species=3
$end

$specie
	name='ions-base-He{He}'
	T=1,
	m=2.0
	q=2.0
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'DF/fu_iz_TAU{tau}.dat'
	injection_method = 0
	injection_rate_relative = .true.
	injection_rate = {round(int(He)/100,2)}
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end

$specie
	name='ions-base-H{H}'
	T=1,
	m=0.5
	q=1.0
	w = 1.0
	mpi_rank = 0
	motion_method = 0
	mono=.false.
	Umono = 0.0
	injection_file = 'DF/fu_iz_TAU{tau}.dat'
	injection_method = 0
	injection_rate_relative = .true.
	injection_rate = {round(int(H)/100, 2)}
	flow_diag = .true.
	track_positions = .false.
	track_velocity = .false.
	erode_prob = 0.0
	se_prob = 0.0
$end

$specie
	name='electrons-base'
	T=1.0
	m=0.005000000
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
$num_diag_regions
	no_diag_reg = 1
$end

$diag_reg
	y_low=14.0
	y_high=16.0
	z_low=13
	z_high=15
	diag_name = 'pot'
	record_property  = 1
	record_only_steady_state= .false.
	nbins = 1
	record_interval= 1.0
	diag_specie = 1
$end

"""

    with open(f"{name}/{type}T{tau}H{H}He{He}P{pot}.inp", "w+") as output:
        if type == "1":
            output.write(inp_file01)
        else:
            output.write(inp_file02)

taus = ["05", "10", "20", "35", "50"]
Hs = ["02", "05", "10", "40", "60", "90", "95", "98"]
Hes = ["98", "95", "90", "60", "40", "10", "05", "02"]
pots = ["01", "03", "06", "10"]
for tau in taus:
    for pot in pots:
        for i in range(len(Hs)):
            write_new_file("1", tau, Hs[i], Hes[i], pot)
            write_new_file("2", tau, Hs[i], Hes[i], pot)