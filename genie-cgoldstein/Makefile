# 12/2020 clean up...such a mess
# changes made in git branch work
# changes made in git branch master

.SUFFIXES: .o .f .f90

NETCDFLIB = /soft/netcdf/4.0.1.intel/lib
NETCDFINC = /soft/netcdf/4.0.1.intel/include
LDFLAGS = -L$(NETCDFLIB) -I$(NETCDFINC) -lnetcdf

.F.o:
	$(FC) $(FLAGS) $(DOPTS) -c $*.F $(LDFLAGS)
.f.o:
	$(FC) $(FLAGS) -c $*.f $(LDFLAGS)
.f:
	$(FC) $(FLAGS) -o $@ $<
.f90.o:
	$(FC) $(FLAGS) -c $*.f90 $(LDFLAGS)
.f90:
	$(FC) $(FLAGS) -o $@ $<

CG  = ../genie-cgoldstein/
BGM = ../genie-biogem/
SGM = ../genie-sedgem/
GLT = ../genie-gemlite/
LND = ../genie-land/
SLD = ../genie-simpleland/
ACM = ../genie-atchem/
GMN = ../genie-main/

# ******************************************
# Object: C-GOLDSTEIN (CG)
OBJS =	$(CG)co.o $(CG)diag.o $(CG)diag2.o $(CG)diagend.o $(CG)diagopsi.o \
	$(CG)diagosc.o $(CG)drgset.o $(CG)inm.o $(CG)invert.o $(CG)island.o \
	$(CG)mains.o $(CG)outm.o $(CG)jbar.o $(CG)ubarsolv.o $(CG)velc.o \
	$(CG)wind.o $(CG)surflux.o $(CG)diaga.o $(CG)diag3.o $(CG)gseto.o \
	$(CG)gseta.o $(CG)tstepa.o $(CG)tstepo.o $(CG)tstipa.o $(CG)tstipo.o \
	$(CG)tstepsic.o $(CG)readroff.o $(CG)radfor.o

# Object: CG + MARINE BIOGEOCHEMISTRY
OBJS_CB = $(OBJS) $(GMN)gem_cmn.o $(GMN)gem_util.o $(GMN)gem_carbchem.o \
	$(ACM)atchem_lib.o $(ACM)atchem_box.o \
	$(ACM)atchem_data.o $(ACM)atchem_main.o \
	$(GMN)cpl_flux_atchem.o $(GMN)cpl_comp_atchem.o \
        $(GMN)genie_netcdf.o \
	$(BGM)biogem_lib.o $(BGM)biogem_box.o \
	$(BGM)biogem_data_netCDF.o $(BGM)biogem_data.o $(BGM)biogem_main.o

# Object: CB + MARINE BIOGEOCHEMISTRY + ENTS
OBJS_CBE = $(OBJS_CB) $(SLD)carbon.o $(SLD)setup_ents.o \
        $(SLD)in_ents.o $(SLD)out_ents.o $(SLD)carbt_diags.o \
        $(SLD)annav_diags.o $(SLD)physt_diags.o $(SLD)screen_diags.o \
        $(SLD)field_interp.o $(SLD)greenland_melt.o $(SLD)ocean_alb.o \
        $(SLD)sealevel.o $(SLD)setup_emissions.o

# Object: CB + MARINE BIOGEOCHEMISTRY + ENTS + SEDGEM
OBJS_CBS = $(OBJS_CB) $(SGM)sedgem_lib.o $(SGM)sedgem_box.o \
	$(SGM)sedgem_data_netCDF.o $(SGM)sedgem_data.o $(SGM)sedgem_main.o \
	$(GMN)cpl_flux_sedgem.o $(GMN)cpl_comp_sedgem.o

# ******************************************
# Source: CB
FBJS =	$(CG)co.f $(CG)diag.f $(CG)diag2.f $(CG)diagend.F $(CG)diagopsi.f \
	$(CG)diagosc.F $(CG)drgset.f $(CG)inm.f $(CG)invert.f $(CG)island.f \
	$(CG)mains.F $(CG)outm.f $(CG)jbar.f $(CG)ubarsolv.f $(CG)velc.f \
	$(CG)wind.f $(CG)surflux.F $(CG)diaga.f $(CG)diag3.f $(CG)gseto.F \
	$(CG)gseta.F $(CG)tstepa.f $(CG)tstepo.F $(CG)tstipa.f $(CG)tstipo.F \
	$(CG)tstepsic.f $(CG)readroff.f $(CG)radfor.F

# Source: CB + MARINE BIOGEOCHEMISTRY
FBJS_CB = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main.f90

# Source: CB + MARINE BIOGEOCHEMISTRY + ENTS
FBJS_CBE = $(FBJS_CB) $(SLD)carbon.F $(SLD)setup_ents.F \
        $(SLD)in_ents.f $(SLD)out_ents.f $(SLD)carbt_diags.f \
        $(SLD)annav_diags.F $(SLD)physt_diags.f $(SLD)screen_diags.f \
        $(SLD)field_interp.f $(SLD)greenland_melt.F $(SLD)ocean_alb.f \
        $(SLD)sealevel.F $(SLD)setup_emissions.f

# Source: CB + MARINE BIOGEOCHEMISTRY + ENTS + SEDGEM
FBJS_CBS = $(FBJS_CB) $(SGM)sedgem_lib.f90  $(SGM)sedgem_box.f90 \
	$(SGM)sedgem_data_netCDF.f90 $(SGM)sedgem_data.f90 $(SGM)sedgem_main.f90 \
	$(GMN)cpl_flux_sedgem.f90 $(GMN)cpl_comp_sedgem.f90

# ******************************************
#### fortran compiler command name
#FC = pgf95
FC = ifort

#### f90 COMPILER FLAGS
# use -g to compile with debugger - otherwise use 2nd row)
FLAGS = -r8 -unroll -O3
#FLAGS = -r8 -unroll -O3 -g

## DOPTS = -D55Ma (55Ma continental config) -D251Ma (251Ma continental config)
DOPTS = -Ddimpa -Ddiso
DOPTS_CB = $(DOPTS) -Dbiogem -Datchem -DSnorm

# wor16b is the original with ez0=0.1; wor16_2in75 uses grid from A. Yool
# MM indicates monthly mean winds (ecmwf)
DOPTS_C16 = $(DOPTS) -Dwor16b
DOPTS_C16_2IN100MM = $(DOPTS) -Dwor16_2in100ecmm

DOPTS_CB16 = $(DOPTS_CB) -Dwor16b
DOPTS_CB16_2IN100 = $(DOPTS_CB) -Dwor16_2in100
DOPTS_CB16_2IN100_neg5deg = $(DOPTS_CB) -Dwor16_2in100 -Dtneg5remin -Dtneg5prod -Dtneg5sol
DOPTS_CB16_2IN100D = $(DOPTS_CB) -Dwor16_2in100d
DOPTS_CB16_2IN100EC = $(DOPTS_CB) -Dwor16_2in100ec
DOPTS_CB16_2IN100MM = $(DOPTS_CB) -Dwor16_2in100ecmm
DOPTS_CB16_2IN100_NODISO = -Ddimpa -Dbiogem -Datchem -DSnorm -Dwor16_2in100

DOPTS_CBE16_2IN100   = $(DOPTS_CB16_2IN100) -Dents
DOPTS_CBE16_2IN100MM   = $(DOPTS_CB16_2IN100MM) -Dents

DOPTS_CBS   = $(DOPTS_CB) -Dsedgem
DOPTS_CBST  = $(DOPTS_CBS) -Dgemlite
DOPTS_CL    = $(DOPTS) -Ddland
DOPTS_CBL   = $(DOPTS_CB) -Ddland
DOPTS_CBSL  = $(DOPTS_CBS) -Ddland
DOPTS_CBSTL = $(DOPTS_CBST) -Ddland

DOPTS_CBS16 = $(DOPTS_CB16) -Dsedgem
#kst 8/18/08:   the following line for 'some reason' doesn't work:  define sedgem in command line for cbs_mesmo  1/28/09:  use caps!
#DOPTS_CBS16_2in100 = $(DOPTS_CB16_2in100) -Dsedgem
DOPTS_CBS16_2IN100 = $(DOPTS_CB16_2IN100) -Dsedgem

CMN = var.cmn Makefile

# ******************************************
# MESMO 1 (cb_2in100_season_ssmax1_kv shortened to cb_2in10_season_kv)

# c goldstein (ocean + embm)
# sent to Andy Price for tuning: c_2in100_seaswind:  the physics only version
c_2in100_seaswind: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_C16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o c_2in100_seaswind $(LDFLAGS)

# basic cb (ocean + embm + biogem)
# definitions of ifdefs
#    dosc = seasonal insolation
#    atmdiff = changed atmospheric thermal diffusivity
#    river = riverine input
#    constBio = constant Bio, used to initialize
#    kvprof = uses a tanh profile of kv rather than a constant
#    DOPTS_CB16_2in75D : uses the dredged Drake passage grid, all others defined in command line
#    nodiso = no isoneutral diffusivity, else sHHHHHHHHHHHHHsmax = 1,5,10
#    ssmax1 = GM parameterization coefficient reduction to keep the model from blowing up - see MESMO1 paper
#    
# SMESMO = MESMO with seasonal ecmwf winds:     dopts_cb16_2in100mm = monthly wind field
# cb_smesmo = cb_2in100_seaswind + utusurf = cb_2in100_seaswind_kv+utusurf
#   uses ecmwf analyzed surface winds (utvel) for usurf=>gasex(ws) this requires the -Dutusurf switch

cb_mesmo: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CB) $(CG)kvprofile.f90 -o cb_mesmo $(LDFLAGS)

cb_smesmo: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo $(LDFLAGS)

cb_smesmo_powsi2n: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dpowsi2n $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_powsi2n $(LDFLAGS)

cbe_mesmo: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo $(LDFLAGS)

mesmo1_ents: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dpowsi2n -Dalb965 -Dcisotopes_ents $(FBJS_CBE) $(CG)kvprofile.f90 -o mesmo1_ents $(LDFLAGS)

cbe_smesmo: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cbe_smesmo $(LDFLAGS)

cbs_mesmo: $(FBJS_CBS) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Dsedgem -Ddosc -Dssmax1 -Dkvprof $(FBJS_CBS) $(CG)kvprofile.f90 -o cbs_mesmo $(LDFLAGS)

# ******************************************
# KM - Nov, 2011 - create four models: MESMO2, MESMOG, MESMO2-ENTS, and MESMOG-ENTS
# have to always have -Dlndice to have ice mask
#
# mesmo2 (mesmo w/ monthly ecmwf winds, icemask PD, ocean biogeo w/ iron, silica, diatoms)
# 	it is essentially cb_smesmo_PDice_tausurf_powsi2n
#
# MESMO2-ENTS
#	add alb965 (albedo reduced by 3.5% to match SST), add C isooptes on land
#	make sure to designate icemask.dat (same as icePD) in goin_bioents
#	-Dewinds overrides ENTS ("new NCEP") winds with MESMO2 (ECMWF) winds; both seasonal but facilitates comparison with MESMO2 
#	should be: cbe_smesmo_alb965_powsi2n_entsCisotopes_sDewinds
#
# MESMOG (mesmo w/ glacial BC)
# 	-Dlgm desinates both icemask and orbital config for lgm
#	-Dicewind -Dice100 designate lgm winds
#
# 	something like the following:
#	mesmog: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
#	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dpowsi2n -Dlndice -Dlgm -Dicewind -Dice100 $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmog $(LDFLAGS)
#
# MESMOG-ENTS
#	add the same alb965 factor
#	make sure to designate icemask_lgm.dat (same as iceLGM) in goin_bioentsa
#	-Dewinds overrides ENTS winds (esp. now that MESMOG winds have glacial anomalies added)
#
# 	something like the following:
#	mesmog_ents: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
#	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dpowsi2n -Dalb965 -Dcisotopes_ents -Dewinds -Dlgm -Dicewind -Dice100 $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmog_ents $(LDFLAGS)
#
#
# KM - June 2012 - finally created MSEMO2 and MESMO2E (see GMD submission)
# From here on out, mesmo2 and mesmo2_ents refer to the following:
# MESMO2  = mesmo2_alb965_ssmaxZ2_nopowsi2n
# MESMO2_ENTS = mesmo2_ents_alb945_ssmaxZ2_nopowsi2n
#
# Some notes on ifdefs:
#    albXXX = reduce global albedo to XXX (e.g., alb965 means 3.5% reduction; see surflux.F); this is taken out of ents ifdef, so works even w/o ents
#    swinds = scaling of wind stress in Southern Hemisphere; need to change that scaling hardwired in BOTH gseta.F and gseto.F
#    ssmaxZX = Gent-McWilliams parameterization coeffcient is depth dependent instead of being constant; see email from Kevin Oliver
#    powsi2n = si2n power law dependence on FeT; not needed in MEMSO2 since it is too weak; inverse relation is hardwired in biogem_box.m
#    cisotopes_ents = adds C13 and C14 in ents; obviously only works when ents is coupled
#    tneg5sol = 5C cooling applied only to gas solubility 
#    fe_so = activates sprinkle_fe, a scaling to increase Fe input in SO; need to change scaling hardwired in biogem_data.f90
#
mesmo2: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2 $(LDFLAGS)

mesmo2_ents: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dcisotopes_ents -Dewinds $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents $(LDFLAGS)


# ******************************************
# KM - April, 2021
# MESMO 3 (includes DOMr, which can be activated in gem_xxx_config.par

mesmo3: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o mesmo3 $(LDFLAGS)

m3_lgmice: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -Dlgmice -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgmice $(LDFLAGS)

m3_lgmorbit: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dorbitlgm -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgmorbit $(LDFLAGS)

m3_lgmrad: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dlgmrad -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgmrad $(LDFLAGS)

# for full lgm bc, need also SO winds, which is controlled by a goin parameter
m3_lgm: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -Dlgmice -Dorbitlgm -Dlgmrad -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgm $(LDFLAGS)
# ***************************************

clean:
	rm -f *.o *.mod

# dependencies

$(CG)co.o: $(CG)co.f $(CMN)
$(CG)diag.o: $(CG)diag.f $(CMN)
$(CG)diag2.o: $(CG)diag2.f $(CMN)
$(CG)diag3.o: $(CG)diag3.f $(CMN)
$(CG)diaga.o: $(CG)diaga.f $(CMN)
$(CG)diagend.o: $(CG)diagend.F $(CMN)
$(CG)diagopsi.o: $(CG)diagopsi.f $(CMN)
$(CG)diagosc.o: $(CG)diagosc.F $(CMN)
$(CG)drgset.o: $(CG)drgset.f $(CMN)
$(CG)inm.o: $(CG)inm.f $(CMN)
$(CG)invert.o: $(CG)invert.f $(CMN)
$(CG)island.o: $(CG)island.f $(CMN)
$(CG)mains.o: $(CG)mains.F $(CMN)
$(CG)outm.o: $(CG)outm.f $(CMN)
$(CG)jbar.o: $(CG)jbar.f $(CMN)
$(CG)ubarsolv.o: $(CG)ubarsolv.f $(CMN)
$(CG)velc.o: $(CG)velc.f $(CMN)
$(CG)wind.o: $(CG)wind.f $(CMN)
$(CG)surflux.o: $(CG)surflux.F $(CMN)
$(CG)gseto.o: $(CG)gseto.F $(CMN)
$(CG)gseta.o: $(CG)gseta.F $(CMN)
$(CG)tstepa.o: $(CG)tstepa.f $(CMN)
$(CG)tstepo.o: $(CG)tstepo.F $(CMN)
$(CG)tstipa.o: $(CG)tstipa.f $(CMN)
$(CG)tstipo.o: $(CG)tstipo.F $(CMN)
$(CG)tstipsic.o: tstipsic.F $(CMN)
$(CG)readroff.o: $(CG)readroff.f $(CMN)
$(CG)radfor.o: $(CG)radfor.F $(CMN)

