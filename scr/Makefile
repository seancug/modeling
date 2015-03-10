output = ../bin
OPT=-O3 -Wall
mcc=g++

objects = modeling.o \
absorb.o \
add_sour.o \
AEA_planer_update_s.o\
AEA_planer_update_v.o\
AEA_topo_update_s.o\
AEA_topo_update_v.o\
allocateMatrix.o\
caculate_coff.o \
checkFD.o \
check_par.o \
flag_mark.o \
Model.o \
outputdata.o\
outseis.o \
pml_coff.o \
rd_model.o \
rd_sour.o \
read_Par.o \
receiver.o \
saveSeis.o \
seismogram.o \
slope_model.o \
snap.o \
surface.o \
topo_surface.o \
topo_update_s.o \
topo_update_v.o \
update_s.o \
update_v.o \
writedsk.o \
Write_Dedault_Seg_Y_head.o \
write_repos.o

modeling:$(objects)
	$(mcc) $(OPT) -o $(output)/modeling $^
modeling.o:modeling.cpp topoft.h glo_var.h
	$(mcc) $(OPT) -c $<
absorb.o:absorb.cpp topoft.h
	$(mcc) $(OPT) -c $<
add_sour.o : add_sour.cpp topoft.h
	$(mcc) $(OPT) -c $<
AEA_planer_update_s.o:AEA_planer_update_s.cpp topoft.h
	$(mcc) $(OPT) -c $<
AEA_planer_update_v.o:AEA_planer_update_v.cpp topoft.h
	$(mcc) $(OPT) -c $<
AEA_topo_update_s.o:AEA_topo_update_s.cpp topoft.h
	$(mcc) $(OPT) -c $<
AEA_topo_update_v.o:AEA_topo_update_v.cpp topoft.h
	$(mcc) $(OPT) -c $<
allocateMatrix.o:allocateMatrix.cpp topoft.h
	$(mcc) $(OPT) -c $<
caculate_coff.o:caculate_coff.cpp topoft.h
	$(mcc) $(OPT) -c $<
checkFD.o:checkFD.cpp topoft.h
	$(mcc) $(OPT) -c $<
check_par.o:check_par.cpp topoft.h
	$(mcc) $(OPT) -c $<
flag_mark.o:flag_mark.cpp topoft.h
	$(mcc) $(OPT) -c $<
Model.o:Model.cpp topoft.h
	$(mcc) $(OPT) -c $<
outputdata.o:outputdata.cpp topoft.h
	$(mcc) $(OPT) -c $<
outseis.o:outseis.cpp topoft.h segy.h
	$(mcc) $(OPT) -c $<
pml_coff.o:pml_coff.cpp topoft.h
	$(mcc) $(OPT) -c $<
rd_model.o:rd_model.cpp topoft.h
	$(mcc) $(OPT) -c $<
rd_sour.o:rd_sour.cpp topoft.h
	$(mcc) $(OPT) -c $<
read_Par.o:read_Par.cpp topoft.h
	$(mcc) $(OPT) -c $<
receiver.o:receiver.cpp topoft.h
	$(mcc) $(OPT) -c $<
saveSeis.o:saveSeis.cpp topoft.h
	$(mcc) $(OPT) -c $<
seismogram.o:seismogram.cpp topoft.h
	$(mcc) $(OPT) -c $<
slope_model.o:slope_model.cpp topoft.h
	$(mcc) $(OPT) -c $<
snap.o:snap.cpp topoft.h
	$(mcc) $(OPT) -c $<
surface.o:surface.cpp topoft.h
	$(mcc) $(OPT) -c $<
topo_surface.o:topo_surface.cpp topoft.h
	$(mcc) $(OPT) -c $<
topo_update_s.o:topo_update_s.cpp topoft.h
	$(mcc) $(OPT) -c $< 
topo_update_v.o:topo_update_v.cpp topoft.h
	$(mcc) $(OPT) -c $< 
update_s.o:update_s.cpp topoft.h
	$(mcc) $(OPT) -c $<
update_v.o:update_v.cpp topoft.h
	$(mcc) $(OPT) -c $< 
writedsk.o:writedsk.cpp topoft.h
	$(mcc) $(OPT) -c $< 
Write_Dedault_Seg_Y_head.o:Write_Dedault_Seg_Y_head.cpp topoft.h
	$(mcc) $(OPT) -c $< 
write_repos.o:write_repos.cpp topoft.h
	$(mcc) $(OPT) -c $<

clean:
	del $(objects) modeling *~ *.gch
