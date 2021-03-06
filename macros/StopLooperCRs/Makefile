CC = g++
ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
INCLUDE = -I../ -I./ $(ROOFITINCLUDE)
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags) $(shell root-config --libs) -lMinuit -lMathMore
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector -lMathMore
endif

SOURCES = ../Core/PartonCombinatorics.cc StopTreeLooper.cc ../Plotting/PlotUtilities.cc ../Core/mt2w_bisect.cc ../Core/mt2bl_bisect.cc ../Core/MT2Utility.cc ../Core/stopUtils.cc ../../CORE/Thrust.cc ../../CORE/EventShape.cc ../Core/STOPT.cc ../../Tools/BTagReshaping/BTagReshaping.cc ../../Tools/BTagReshaping/btag_payload_light.cc ../../Tools/BTagReshaping/btag_payload_b.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libStopTreeLooper.so

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@  $(shell root-config --libs) -lMinuit

LinkDef_out.cxx: LinkDef.h StopTreeLooper.h ../Plotting/PlotUtilities.h ../Core/stopUtils.h   
	rootcint -f $@ -c $(INCLUDE) StopTreeLooper.h ../Plotting/PlotUtilities.h ../Core/stopUtils.h $<

# General rule for making object files
%.d:	%.cc
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.d:	%.cxx
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

LIBS = $(LIB) 

.PHONY: all
all:	$(LIBS)  

.PHONY: clean
clean:  
	rm -f *.d \
	rm -f *.o \
	rm -f */*.d \
	rm -f */*.o \
	rm -f *.so \
	rm -f LinkDef_out* 

-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

