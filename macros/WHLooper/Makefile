CC = g++
ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#INCLUDE = -I../ -I./ $(ROOFITINCLUDE)
INCLUDE = -I../ -I./ 
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags) $(shell root-config --libs) -lMathMore
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector -lMathMore
endif

SOURCES = WHLooper.cc ../../Tools/BTagReshaping/BTagReshaping.cc ../../Tools/BTagReshaping/btag_payload_light.cc ../../Tools/BTagReshaping/btag_payload_b.cc ../Plotting/PlotUtilities.cc ../Core/PartonCombinatorics.cc ../Core/MT2Utility.cc ../Core/mt2w_bisect.cc ../Core/mt2bl_bisect.cc ../Core/MT2.cc ../Core/stopUtils.cc ../Core/STOPT.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libWHLooper.so

CORESOURCES = ../../CORE/CMS2.cc ../../CORE/utilities.cc ../../CORE/ssSelections.cc ../../CORE/electronSelections.cc ../../CORE/electronSelectionsParameters.cc ../../CORE/MITConversionUtilities.cc ../../CORE/muonSelections.cc ../../CORE/eventSelections.cc ../../CORE/trackSelections.cc ../../CORE/metSelections.cc ../../CORE/jetSelections.cc ../../CORE/photonSelections.cc ../../CORE/triggerUtils.cc ../../CORE/triggerSuperModel.cc ../../CORE/mcSelections.cc ../../CORE/susySelections.cc ../../CORE/mcSUSYkfactor.cc ../../CORE/SimpleFakeRate.cc ../../Tools/goodrun.cc ../../Tools/vtxreweight.cc ../../Tools/msugraCrossSection.cc  ../../CORE/jetsmear/JetSmearer.cc ../../CORE/jetsmear/JetResolution.cc ../../CORE/jetsmear/SigInputObj.cc ../../CORE/jetSmearingTools.cc ../../CORE/QuarkGluonTagger/QGLikelihoodCalculator.cc ../../CORE/QuarkGluonTagger/QuarkGluonTagger.cc 

COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libWHLooperCORE.so

LIBS = $(LIB) $(CORELIB)

libs:	$(LIBS)

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 

$(CORELIB):	$(COREOBJECTS) 
	echo "Linking $(CORELIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

LinkDef_out.cxx: LinkDef.h WHLooper.h ../Plotting/PlotUtilities.h ../Core/stopUtils.h 
	rootcint -f $@ -c $(INCLUDE) WHLooper.h ../Plotting/PlotUtilities.h ../Core/stopUtils.h $<

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
	rm -f ../Core/*.d \
	rm -f ../Core/*.o \
	rm -f LinkDef_out* 

-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

