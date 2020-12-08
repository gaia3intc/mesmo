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

OBJS =	$(CG)co.o $(CG)diag.o $(CG)diag2.o $(CG)diagend.o $(CG)diagopsi.o \
	$(CG)diagosc.o $(CG)drgset.o $(CG)inm.o $(CG)invert.o $(CG)island.o \
	$(CG)mains.o $(CG)outm.o $(CG)jbar.o $(CG)ubarsolv.o $(CG)velc.o \
	$(CG)wind.o $(CG)surflux.o $(CG)diaga.o $(CG)diag3.o $(CG)gseto.o \
	$(CG)gseta.o $(CG)tstepa.o $(CG)tstepo.o $(CG)tstipa.o $(CG)tstipo.o \
	$(CG)tstepsic.o $(CG)readroff.o $(CG)radfor.o

OBJS_CB = $(OBJS) $(GMN)gem_cmn.o $(GMN)gem_util.o $(GMN)gem_carbchem.o \
	$(ACM)atchem_lib.o $(ACM)atchem_box.o \
	$(ACM)atchem_data.o $(ACM)atchem_main.o \
	$(GMN)cpl_flux_atchem.o $(GMN)cpl_comp_atchem.o \
        $(GMN)genie_netcdf.o \
	$(BGM)biogem_lib.o $(BGM)biogem_box.o \
	$(BGM)biogem_data_netCDF.o $(BGM)biogem_data.o $(BGM)biogem_main.o

OBJS_CBE = $(OBJS_CB) $(SLD)carbon.o $(SLD)setup_ents.o \
        $(SLD)in_ents.o $(SLD)out_ents.o $(SLD)carbt_diags.o \
        $(SLD)annav_diags.o $(SLD)physt_diags.o $(SLD)screen_diags.o \
        $(SLD)field_interp.o $(SLD)greenland_melt.o $(SLD)ocean_alb.o \
        $(SLD)sealevel.o $(SLD)setup_emissions.o

OBJS_CBS = $(OBJS_CB) $(SGM)sedgem_lib.o $(SGM)sedgem_box.o \
	$(SGM)sedgem_data_netCDF.o $(SGM)sedgem_data.o $(SGM)sedgem_main.o \
	$(GMN)cpl_flux_sedgem.o $(GMN)cpl_comp_sedgem.o

OBJS_CBST = $(OBJS_CBS) $(GLT)gemlite_lib.o $(GLT)gemlite_box.o \
	$(GLT)gemlite_data.o $(GLT)gemlite_main.o

OBJS_L = $(LND)land_coupler.o $(LND)canopy.o $(LND)growth.o \
	$(LND)leaf_c4.o $(LND)lotka.o $(LND)pft_sparm.o $(LND)raero.o \
	$(LND)swrad.o $(LND)compete.o $(LND)land.o $(LND)leaf_c3.o \
	$(LND)microbe.o $(LND)phenol.o $(LND)sf_stom.o $(LND)triffid.o \
	$(LND)decay.o $(LND)land_state.o $(LND)leaf_lit.o $(LND)penmon.o \
	$(LND)qsat.o $(LND)soilcarb.o $(LND)vegcarb.o $(LND)inml.o \
	$(LND)outml.o $(LND)diagend_land.o $(LND)gsetl.o $(LND)diag_land.o

OBJS_CBL = $(OBJS_CB) $(OBJS_L)
OBJS_CBSTL = $(OBJS_CBST) $(OBJS_L)
OBJS_CBSL = $(OBJS_CBS) $(OBJS_L)
OBJS_CL = $(OBJS) $(OBJS_L)

FBJS =	$(CG)co.f $(CG)diag.f $(CG)diag2.f $(CG)diagend.F $(CG)diagopsi.f \
	$(CG)diagosc.F $(CG)drgset.f $(CG)inm.f $(CG)invert.f $(CG)island.f \
	$(CG)mains.F $(CG)outm.f $(CG)jbar.f $(CG)ubarsolv.f $(CG)velc.f \
	$(CG)wind.f $(CG)surflux.F $(CG)diaga.f $(CG)diag3.f $(CG)gseto.F \
	$(CG)gseta.F $(CG)tstepa.f $(CG)tstepo.F $(CG)tstipa.f $(CG)tstipo.F \
	$(CG)tstepsic.f $(CG)readroff.f $(CG)radfor.F

##moved biogem_lib.f90 and biogem_box.f90 up to line 2, so they'd be available to atchm_main.90...     ORG<KEEP>
FBJS_CB = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main.f90
## FBJS_CB_xxx to test with different biogem f90files TaTa 150929

##FBJS_CB_2N1T for 2N1T simple model
FBJS_CB_2N1T = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib_2N1T.f90 $(BGM)biogem_box_2N1T.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF_2N1T.f90 $(BGM)biogem_data_2N1T.f90 $(BGM)biogem_main_2N1T.f90

## FBJS_CB_jz changes for JZ's DOM stuff Tata 180326
FBJS_CB_jz = $(FBJS) $(GMN)gem_cmn_jz.f90 $(GMN)gem_util_jz.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib_jz.f90 $(BGM)biogem_box_jz.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF_jz.f90 $(BGM)biogem_data_jz.f90 $(BGM)biogem_main_jz.f90

FBJS_CB_c = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main_c.f90

FBJS_CB_c2 = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main_c2.f90

FBJS_CB_d = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main_d.f90

FBJS_CB_e = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main_e.f90

FBJS_CBTMP = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box_noTreminDOM.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main.f90

FBJS_CBRS = $(FBJS) $(GMN)gem_cmn.f90 $(GMN)gem_util.f90 $(GMN)gem_carbchem.f90 \
	$(ACM)atchem_lib.f90 $(ACM)atchem_box.f90 \
	$(ACM)atchem_data.f90 $(ACM)atchem_main.f90 \
	$(GMN)cpl_flux_atchem.f90 $(GMN)cpl_comp_atchem.f90 \
        $(GMN)genie_netcdf.f90 \
	$(BGM)biogem_lib.f90 $(BGM)biogem_box2Nrestore.f90 \
	$(BGM)biogem_data_netCDF.f90 $(BGM)biogem_data.f90 $(BGM)biogem_main.f90

FBJS_CBE = $(FBJS_CB) $(SLD)carbon.F $(SLD)setup_ents.F \
        $(SLD)in_ents.f $(SLD)out_ents.f $(SLD)carbt_diags.f \
        $(SLD)annav_diags.F $(SLD)physt_diags.f $(SLD)screen_diags.f \
        $(SLD)field_interp.f $(SLD)greenland_melt.F $(SLD)ocean_alb.f \
        $(SLD)sealevel.F $(SLD)setup_emissions.f

FBJS_CBS = $(FBJS_CB) $(SGM)sedgem_lib.f90  $(SGM)sedgem_box.f90 \
	$(SGM)sedgem_data_netCDF.f90 $(SGM)sedgem_data.f90 $(SGM)sedgem_main.f90 \
	$(GMN)cpl_flux_sedgem.f90 $(GMN)cpl_comp_sedgem.f90

FBJS_CBST = $(FBJS_CBS) $(GLT)gemlite_lib.f90 $(GLT)gemlite_box.f90 \
	$(GLT)gemlite_data.f90 $(GLT)gemlite_main.f90

FBJS_L = $(LND)land_coupler.F $(LND)canopy.f $(LND)growth.f \
	$(LND)leaf_c4.f $(LND)lotka.f $(LND)pft_sparm.f $(LND)raero.f \
	$(LND)swrad.f $(LND)compete.f $(LND)land.f $(LND)leaf_c3.f \
	$(LND)microbe.f $(LND)phenol.f $(LND)sf_stom.f $(LND)triffid.f \
	$(LND)decay.f $(LND)land_state.f $(LND)leaf_lit.f $(LND)penmon.f \
	$(LND)qsat.f $(LND)soilcarb.f $(LND)vegcarb.f $(LND)inml.f \
	$(LND)outml.f $(LND)diagend_land.f $(LND)gsetl.f $(LND)diag_land.f

FBJS_CBL = $(FBJS_CB) $(FBJS_L)
FBJS_CBSTL = $(FBJS_CBST) $(FBJS_L)
FBJS_CBSL = $(FBJS_CBS) $(FBJS_L) 
FBJS_CL = $(FBJS) $(FBJS_L)

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

# by M.Chikamoto 07-24-2006
DOPTS_CBN = $(DOPTS_CB) -DfixN
DOPTS_CBS_N = $(DOPTS_CBS) -DfixN
DOPTS_CBS_N14C = $(DOPTS_CBS) -DfixN -Dprod14C

CMN = var.cmn Makefile

# c goldstein (ocean + embm)
#----sent to Andy Price for tuning: c_2in100_seaswind:  the physics only version
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

# --------------- 2in100 models:----------------------------------------------

#-----tuning variations:
cb_2in100_dredge: $(FBJS_CB) $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100D) $(FBJS_CB) -o cb_2in100_dredge $(LDFLAGS)
cb_2in100_dredge_constBio: $(FBJS_CB) $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100D) -DconstBio $(FBJS_CB) -o cb_2in100_dredge_constBio $(LDFLAGS)

#%%%%%%%%--------MESMO1             (cb_2in100_season_ssmax1_kv shortened to cb_2in10_season_kv:) shortened to cb_mesmo  THIS IS MESMO  ***************
cb_2in100_season_ssmax1_kv_constBio: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Ddosc -Dssmax1 -Dkvprof -DconstBio $(FBJS_CB) $(CG)kvprofile.f90 -o cb_2in100_season_ssmax1_kv_constBio $(LDFLAGS)

cb_mesmo: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CB) $(CG)kvprofile.f90 -o cb_mesmo $(LDFLAGS)

cb_mesmo_lgm: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dlgm -Dlndice $(FBJS_CB) $(CG)kvprofile.f90 -o cb_mesmo_lgm $(LDFLAGS)

cb_mesmo_powsi2n_PDice: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dpowsi2n -Dlndice -DicePD $(FBJS_CB) $(CG)kvprofile.f90 -o cb_mesmo_powsi2n_PDice $(LDFLAGS)

cb_mesmo_ema: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Demaice $(FBJS_CB) $(CG)kvprofile.f90 -o cb_mesmo_ema $(LDFLAGS)

#%%%%%%%%%%%%  MESMO1 with ENTS:
cbe_mesmo: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo $(LDFLAGS)

cbe_mesmo_cldalb: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dcldalb $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_cldalb $(LDFLAGS)

cbe_mesmo_alb965_powsi2n_entsCisotopes: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dpowsi2n -Dcisotopes_ents $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_alb965_powsi2n_entsCisotopes $(LDFLAGS)

mesmo1_ents: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dpowsi2n -Dalb965 -Dcisotopes_ents $(FBJS_CBE) $(CG)kvprofile.f90 -o mesmo1_ents $(LDFLAGS)

mesmo1_ents_ppm2gtc_joos: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dpowsi2n -Dalb965 -Dcisotopes_ents $(FBJS_CBE) $(CG)kvprofile.f90 -o mesmo1_ents_ppm2gtc_joos $(LDFLAGS)

cbe_mesmo_alb965_powsi2n_entsCisotopes_tneg5sol: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dpowsi2n -Dcisotopes_ents -Dtnegsol $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_alb965_powsi2n_entsCisotopes_tneg5sol $(LDFLAGS)

cbe_mesmo_alb965_powsi2n_entsCisotopes_feso: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dpowsi2n -Dcisotopes_ents -Dfe_so $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_alb965_powsi2n_entsCisotopes_feso $(LDFLAGS)

cbe_smesmo_alb965_powsi2n_entsCisotopes_sDewinds: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dpowsi2n -Dalb965 -Dcisotopes_ents -Dewinds $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cbe_smesmo_alb965_powsi2n_entsCisotopes_sDewinds $(LDFLAGS)

cbe_mesmo_alb965_powsi2n_entsCisotopes_lgm: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dpowsi2n -Dcisotopes_ents -Dlgm $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_alb965_powsi2n_entsCisotopes_lgm $(LDFLAGS)

cbe_mesmo_alb965_powsi2n_entsCisotopes_10k: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dpowsi2n -Dcisotopes_ents -Dlgm10k $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_alb965_powsi2n_entsCisotopes_10k $(LDFLAGS)

cbe_mesmo_EMIC_init: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dpowsi2n -Demic5_850 $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmo_EMIC_init $(LDFLAGS)

cbe_mesmox: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmox $(LDFLAGS)

cbe_mesmoalbo: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalbox $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmoalbo $(LDFLAGS)

cbe_mesmoalbs: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100) -Ddosc -Dssmax1 -Dkvprof -Dalbsx $(FBJS_CBE) $(CG)kvprofile.f90 -o cbe_mesmoalbs $(LDFLAGS)

cbe_smesmo: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cbe_smesmo $(LDFLAGS)

cbe_smesmo_alb965: $(FBJS_CBE) $(CG)kvprofile.f90 $(CMN) $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(CBE16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dalb965 -Dutusurf $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cbe_smesmo_alb965 $(LDFLAGS)

#%%%%%%%%---------MESMO with sediments:   NOTE:  cbs_2in100_season_kv shortened to cbs_mesmo
cbs_mesmo: $(FBJS_CBS) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Dsedgem -Ddosc -Dssmax1 -Dkvprof $(FBJS_CBS) $(CG)kvprofile.f90 -o cbs_mesmo $(LDFLAGS)

cbs_mesmo_constBio: $(FBJS_CBS) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100) -Dsedgem -Ddosc -Dssmax1 -Dkvprof -DconstBio $(FBJS_CBS) $(CG)kvprofile.f90 -o cbs_mesmo_constBio $(LDFLAGS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

#%%%%%%-----------MESMO with ecmwf annual winds ---_ec & dopts_cb16_2in100ec = mean annual ecmwf 
cb_2in100_seas_kv_ec_constBio: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100EC) -Ddosc -Dssmax1 -Dkvprof -DconstBio $(FBJS_CB) $(CG)kvprofile.f90 -o cb_2in100_seas_kv_ec_constBio $(LDFLAGS)
cb_2in100_seas_kv_ec: $(FBJS_CB) $(CG)kvprofile.f90 $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100EC) -Ddosc -Dssmax1 -Dkvprof $(FBJS_CB) $(CG)kvprofile.f90 -o cb_2in100_seas_kv_ec $(LDFLAGS)

#%%%%%%%%---------SMESMO = MESMO with seasonal ecmwf winds:     dopts_cb16_2in100mm = monthly wind field
# cb_smesmo = cb_2in100_seaswind + utusurf = cb_2in100_seaswind_kv+utusurf
#   uses ecmwf analyzed surface winds (utvel) for usurf=>gasex(ws) this requires the -Dutusurf switch

cb_smesmo: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo $(LDFLAGS)

cb_smesmo_powsi2n: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dpowsi2n $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_powsi2n $(LDFLAGS)

###  holds CO3 const for "constant" CaCO3 production###use this for 'constant' production..(pseudo abio)
cb_smesmo_prodcaco3: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -DCO3const -Dprodtotalconst $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_prodcaco3 $(LDFLAGS)
###note:  this skips then entire bio production scheme in box
cb_smesmo_constBio: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -DconstBio -Dssmax1 -Dkvprof -Dutusurf $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_constBio $(LDFLAGS)

cb_smesmo_lgm: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dlndice -Dlgm $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_lgm $(LDFLAGS)

#%%%%%%%%------------SMESMO with some glacial conditions, icewinds, etc.
#                                         -----NOTE: cb_seaswind_icePD == cb_seaswind_PD_ice (8/5/09)
#---                                                  -smesmo_PDice = lndice albedo = 0.7, utvel
#  note:  utusurf should be used.  if icewinds, it will be scaled according to the icewind scaling
#
cb_smesmo_PDice: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dlndice -DicePD $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_PDice $(LDFLAGS)
cb_smesmo_PDice_tausurf: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dlndice -DicePD $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_PDice_tausurf $(LDFLAGS)
#---uses scaled tau proportionality for usurf gasex  
cb_smesmo_anom_50: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dicewind -Dice50 $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_anom_50 $(LDFLAGS)
cb_smesmo_anom_75A: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dicewind -Dice75A $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_anom_75A $(LDFLAGS)

###  holds CO3 const for "constant" CaCO3 production###use this for 'constant' production..(pseudo abio)
cb_smesmo_prodcaco3_anom_75A: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dicewind -Dice75A -DCO3const -Dprodtotalconst $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_prodcaco3_anom_75A $(LDFLAGS)

cb_smesmo_anom_100: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dicewind -Dice100 $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_anom_100 $(LDFLAGS)
cb_smesmo_anom_50_lgm: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dicewind -Dice50 -Dlndice -Dlgm $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_anom_50_lgm $(LDFLAGS)
cb_smesmo_anom_100_lgm: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dutusurf -Dicewind -Dice100 -Dlndice -Dlgm $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_smesmo_anom_100_lgm $(LDFLAGS)


cb_seaswind_anom_100_lgm10k: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Dicewind -Dice100 -Dlndice -Dlgm10k $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_seaswind_anom_100_lgm10k $(LDFLAGS)
cb_2in100_seaswind_ema: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax1 -Dkvprof -Demaice $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f -o cb_2in100_seaswind_ema $(LDFLAGS)


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

#**********************************
mesmo2_feso: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD $(FBJS_CB) -Dalb965 -Dfe_so $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_feso $(LDFLAGS)

mesmo2_swinds: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD $(FBJS_CB) -Dalb965 -Dswinds $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_swinds $(LDFLAGS)

mesmo2_ents: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dcisotopes_ents -Dewinds $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents $(LDFLAGS)
mesmo2_ents_nopac: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dnopac -Dcisotopes_ents -Dewinds $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_nopac $(LDFLAGS)

mesmo2_ents_nout: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -DicePD -Dcisotopes_ents $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_nout $(LDFLAGS)
mesmo2_ents_nout_sc: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -DicePD -Dcisotopes_ents -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_nout_sc $(LDFLAGS)
#**********NOTE: mesmo2_ents_tausc = mesmo2_ents_nout_sc*********************************************8
mesmo2_ents_tausc: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -DicePD -Dcisotopes_ents -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_tausc $(LDFLAGS)
#******************hadly wind contributions to globtair ***************30oct13************************
mesmo2_ents_tausc_lgmice: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -Dlgmice -Dcisotopes_ents -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_tausc_lgmice $(LDFLAGS)

mesmo2_ents_tausc_lgmorbit: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -DicePD -Dlgmorbit -Dcisotopes_ents -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_tausc_lgmorbit $(LDFLAGS)

mesmo2_ents_tausc_lgm: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -Dlgm -Dcisotopes_ents -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_tausc_lgm $(LDFLAGS)

#*****************************
mesmo2_ents_tauscatl: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -DicePD -Dcisotopes_ents -Datlscfonly -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_tauscatl $(LDFLAGS)

mesmo2_ents_ut_sc: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dutrel -Dlndice -DicePD -Dcisotopes_ents -Dsc $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_ut_sc $(LDFLAGS)

mesmo2_ents_nout_nopac: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dlndice -DicePD -Dnopac -Dcisotopes_ents $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_nout_nopac $(LDFLAGS)

mesmo2_ents_hadley: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dhadley -Dlndice -DicePD -Dcisotopes_ents -Dewinds $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_hadley $(LDFLAGS)
#-Dhadly currently wired for 21k hadley winds: (17jul13)
mesmo2_ents_hadley_lgm: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dhadley -Dlndice -Dlgm -Dcisotopes_ents -Dewinds $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_hadley_lgm $(LDFLAGS)

mesmo2_ents_hnzley_lgm: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dhnzley -Dlndice -Dlgm -Dcisotopes_ents -Dewinds $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_hnzley_lgm $(LDFLAGS)

mesmo2_ents_hnzley_noewind_lgm: $(FBJS_CBE) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CBE16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dhnzley -Dlndice -Dlgm -Dcisotopes_ents $(FBJS_CBE) -Dalb945 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_ents_hnzley_noewind_lgm $(LDFLAGS)

# *****************************************
# Tata 180214
# # Testing out gas exchange dumping with mesmo2
# gas_flux_dump = switch to put in gas exchange flux control
mesmo2_gas_dump: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo2_gas_dump $(LDFLAGS)

# ******************************************
# Tata 06/03/15 
# mesmo3 = mesmo2 with simple denitrification and nitrogen fixation scheme
#    notes on ifdefs
#    stoich = Flexible C:N:P:-O2
#    denit_prog = prognostic denitrification
#    gas_flux_dump = switch to put in gas exchange flux control (now default for mesmo3)
#    flex_efratio_henson = flexible export ratio as a function of temperature following Henson11, GRL
#    flex_efratio_laws = flexible export ratio following Laws00, GBC
#
#    Tata 180920 
#    Mesmo3 is the model to use. All the other features can be turned off/on in biogem_config file
#
#    KM 1/2019 std mesmo3 to include Laws flexible ef ratio
#    	W/o Laws or Henson designated, model defaults to Dunne formulation (mesmo3_efd)
#
#    KM 5/2020 m3_docr = previous mesmo3 as described in Paleo paper + DOCr in JZ
#    Jacob's model (control experiment is 901012c)

mesmo3: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o mesmo3 $(LDFLAGS)

m3_docr: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump -Dflex_efratio_laws -Ddocr $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_docr $(LDFLAGS)

m3_nodocr: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_nodocr $(LDFLAGS)

m3_docr_dunne: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump -Ddocr $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_docr_dunne $(LDFLAGS)

m3_dunne: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_dunne $(LDFLAGS)

mesmo3_nogasdamp: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o mesmo3_nogasdamp $(LDFLAGS)
mesmo3_efh: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump -Dflex_efratio_henson $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o mesmo3_efh $(LDFLAGS)
mesmo3_efd: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dstoich -Dgas_flux_dump $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o mesmo3 $(LDFLAGS)

m3_lgmice: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -Dlgmice -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgmice $(LDFLAGS)

m3_lgmorbit: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dorbitlgm -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgmorbit $(LDFLAGS)

m3_lgmrad: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dlgmrad -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgmrad $(LDFLAGS)

# for full lgm bc, need also SO winds, which is controlled by a goin parameter
m3_lgm: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -Dlgmice -Dorbitlgm -Dlgmrad -Dstoich -Dgas_flux_dump -Dflex_efratio_laws $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f $(BGM)stoich.f90 -o m3_lgm $(LDFLAGS)

# *****************************************************************************
# Simplified Version of MESMO (MESMO_simp)
# Tata 180509
# Used for 2N(PO4, NO3) and 1T study
# gas_flux_dump : gas exchange dumping scheme to prevent excessive O2 and CO2 gas exchange
# flex_efratio_dunne = flexible export ratio following Dunne05, GBC
# stoich_2N1T_smodel : power-law flexible C:P and C:N
# stoich_1N1T_smodel : power-law flexible C:P
# stoich_2N1T_GM15 : flexible C:P and C:N of Galbraith and Martiny (2015), PNAS
# stoich_1N1T_GM15 : flexible C:P of Galbraith and Martiny (2015), PNAS
#
mesmo_simp: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump -Dflex_efratio_dunne $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp $(LDFLAGS)
mesmo_simp_stoich: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump -Dflex_efratio_dunne -Dstoich_2N1T_smodel -Dstoich  $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp_stoich $(LDFLAGS)
mesmo_simp_stoich_GM15: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump -Dflex_efratio_dunne -Dstoich_2N1T_GM15 -Dstoich  $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp_stoich_GM15 $(LDFLAGS)
mesmo_simp_stoich_1N1T: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump -Dflex_efratio_dunne -Dstoich_1N1T_smodel -Dstoich  $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp_stoich_1N1T $(LDFLAGS)
mesmo_simp_stoich_1N1T_GM15: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump -Dflex_efratio_dunne -Dstoich_1N1T_GM15 -Dstoich  $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp_stoich_1N1T_GM15 $(LDFLAGS)
mesmo_simp_fixef: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp_fixef $(LDFLAGS)
#mesmo_simp: $(FBJS_CB) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
#	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump $(FBJS_CB) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp $(LDFLAGS)
#mesmo_simp_stoich: $(FBJS_CB_2N1T) $(CG)kvprofile.f90 $(CG)wind_interp.f $(CMN)
#	$(FC) $(FLAGS) $(DOPTS_CB16_2IN100MM) -Ddosc -Dssmax_z2 -Dkvprof -Dutusurf -Dlndice -DicePD -Dgas_flux_dump -Dstoich -Dstoich_2N1T $(FBJS_CB_2N1T) -Dalb965 $(CG)kvprofile.f90 $(CG)wind_interp.f -o mesmo_simp_stoich $(LDFLAGS)


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

