F77 = mpif90 -i4 -real-size 32 -O2

FILES =  dimensions.f global.f grid_interp.f maind.f gutsf.f gutsp.f misc.f boundary.f part_init.f initial.f inputs.f chem_rates.f
INCLUDE = dimensions.f
OBJECTS = dimensions.o inputs.o global.o boundary.o initial.o grid_interp.o gutsf.o misc.o gutsp.o part_init.o chem_rates.o maind.o
MODS = dimensions.mod global.mod inputs.mod boundary.mod initial.mod grid_interp.mod gutsf.mod misc.mod  gutsp.mod  part_init.mod chem_rates.mod

hybrid:	$(OBJECTS) 
	$(F77) -o hybrid $(OBJECTS) 

mods:	$(MODS) 

clean:
	rm *.o hybrid *.out *.mod

dimensions.o:dimensions.f $(INCLUDE);$(F77) -c dimensions.f
global.o:global.f $(INCLUDE);$(F77) -c global.f
maind.o:maind.f $(INCLUDE);$(F77) -c maind.f
gutsf.o:gutsf.f $(INCLUDE);$(F77) -c gutsf.f
gutsp.o:gutsp.f $(INCLUDE);$(F77) -c gutsp.f
misc.o:misc.f $(INCLUDE);$(F77) -c misc.f
boundary.o:boundary.f $(INCLUDE);$(F77) -c boundary.f
part_init.o:part_init.f $(INCLUDE);$(F77) -c part_init.f
initial.o:initial.f $(INCLUDE);$(F77) -c initial.f
inputs.o:inputs.f $(INCLUDE);$(F77) -c inputs.f
grid_interp.o:grid_interp.f $(INCLUDE);$(F77) -c grid_interp.f
chem_rates.o:chem_rates.f $(INCLUDE);$(F77) -c chem_rates.f

dimensions.mod:dimensions.f $(INCLUDE);$(F77) -c dimensions.f
global.mod:global.f $(INCLUDE);$(F77) -c global.f
gutsf.mod:gutsf.f $(INCLUDE);$(F77) -c gutsf.f
gutsp.mod:gutsp.f $(INCLUDE);$(F77) -c gutsp.f
misc.mod:misc.f $(INCLUDE);$(F77) -c misc.f
boundary.mod:boundary.f $(INCLUDE);$(F77) -c boundary.f
part_init.mod:part_init.f $(INCLUDE);$(F77) -c part_init.f
initial.mod:initial.f $(INCLUDE);$(F77) -c initial.f
inputs.mod:inputs.f $(INCLUDE);$(F77) -c inputs.f
grid_interp.mod:grid_interp.f $(INCLUDE);$(F77) -c grid_interp.f
chem_rates.mod:chem_rates.f $(INCLUDE);$(F77) -c chem_rates.f