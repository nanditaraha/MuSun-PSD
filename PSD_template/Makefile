.PHONY: clean cleandisplay



ROOT6 := $(shell command -v rootcling 2> /dev/null)
UIUCDIR = $(HOME)/musun_cvs
MTA = $(HOME)/musun_cvs

UIUC_INCLUDE_DIRS = $(UIUCDIR):$(UIUCDIR)/Display:$(UIUCDIR)/Raw:$(UIUCDIR)/Muon:$(UIUCDIR)/Electron:$(UIUCDIR)/Example:$(UIUCDIR)/Clustering:$(UIUCDIR)/MuonElectron/:$(UIUCDIR)/Lifetime:$(UIUCDIR)/Calibration:$(UIUCDIR)/DisplayTrigger:$(UIUCDIR)/Spark:$(UIUCDIR)/Francoise:$(UIUCDIR)/Efficiencies:$(UIUCDIR)/Background:$(UIUCDIR)/Fadc:$(UIUCDIR)/HotWires:$(UIUCDIR)/EventTree:$(UIUCDIR)/macros/EventTree/src:$(UIUCDIR)/Fadc:$(UIUCDIR)/../uky:$(UIUCDIR)/../shared:$(UIUCDIR)/TPC:$(UIUCDIR)/Parameters:$(UIUCDIR)/Neutron:$(UIUCDIR)/MCTruth

UIUC_CFLAGS = $(patsubst %,-I%,$(subst :, ,$(UIUC_INCLUDE_DIRS)))
UIUC_CFLAGS += -DNOMUPC2
UIUC_CFLAGS += -DUIUCTREEANA
UIUC_CFLAGS += -DEPCSTUDY
#UIUC_CFLAGS += -DFADCRUN8

# Use both of these lines for RUN9+, otherwise comment out both
UIUC_CFLAGS += -DRUN9
UIUC_CFLAGS += -DCDEB
#UIUC_CFLAGS += -DWFD
#BIGCATS = $(HOME)/mu/src/uiuc
#UIUC_CFLAGS += -DKICKTIMING
#UIUC_CLFAGS += -DKICKREQ

LIB_DIR = $(MUTRUNK)/work
#LIB_DIR = $(MU)

# Include the lists of source files
ifdef BIGCATS
include ${MTA}/src/files.make.bigcats
else
include ${MTA}/files.make
endif 

ROOT_CFLAGS = `$(ROOTSYS)/bin/root-config --cflags`
# ROOT_LIBS = `$(ROOTSYS)/bin/root-config --libs` -lThread
ROOT_LIBS = `$(ROOTSYS)/bin/root-config --libs` -lThread -lGui
#$(info $(ROOT_CFLAGS))

# Use "make cleandisplay && make" to recompile for event display use.
# uncomment the next line ONLY if you are running the event display
#USE_EVENT_DISPLAY = 1
ifdef USE_EVENT_DISPLAY
  UIUC_CFLAGS += -DEVENTDISPLAY
  SOURCES += EventDisplayTrigger.C EventDisplaySandbox.C
  HEADERS += EventDisplayTrigger.h EventDisplaySandbox.h
endif

# uncomment the next 3 lines ONLY if you are running the BU event display (all detectors)
# also change the choice of display in main.cpp around line 40 and around line 120
#USE_BUEVENT_DISPLAY = 1
ifdef USE_BUEVENT_DISPLAY
  UIUC_CFLAGS += -DEVENTDISPLAY
  UIUC_CFLAGS += -DEVENTDISPLAYBU
  SOURCES += BUEventViewer.C EventDisplaySandbox.C
  HEADERS += BUEventViewer.h EventDisplaySandbox.h
  ROOT_LIBS += $(ROOTSYS)/lib/libEve.so $(ROOTSYS)/lib/libRGL.so $(ROOTSYS)/lib/libGeom.so $(ROOTSYS)/lib/libGeomPainter.so
endif

# show analysis objects on display (line fit to MuStop, stop point)
# UIUC_CFLAGS += -DANASHAPES
# UIUC_CFLAGS += -DMUDISPAUX

SOURCES += $(DICTIONARY)
LIFETIMEHEADERS = LifetimeEvent.h LifetimeEventTreeBuilder.h TLifetimeHistograms.h TLifetimeHistograms_ElectronBG.h TLifetimeHistograms_Tree.h TTPCDelayedEventHistograms.h TLifetimeHistograms_ePC.h LinkDefLifetimeEvent.h
STOPDISTHEADERS = stopDistPair.h stopDistPairTreeBuilder.h LinkDefStopDistTree.h

all: mta libMuETree.so libMuETree.rootmap libDelayedEvent.so libDelayedEvent.rootmap libLifetimeEvent.so libLifetimeEvent.rootmap libStopDistTree.so libStopDistTree.rootmap

# Determine the names of the object files by pattern substitution
# from the source file names.
OBJECTS := $(notdir $(patsubst %.C,%.o,$(SOURCES)))
OBJECTS := $(notdir $(patsubst %.cpp,%.o,$(OBJECTS)))

vpath % .:../src/ucb:$(UIUC_INCLUDE_DIRS):../src/shared:../src/uky:

ifdef CERN
CERN_LIBS   = -L$(CERN)/$(CERN_LEVEL)/lib -lphtools -lgraflib -lgrafX11 \
                                       -lmathlib -lpacklib -lpawlib -lnss_nis -lg2c
else
CERN_LIBS = 
endif

#X11LIBS  = -L/usr/X11R6/lib -lX11  
X11LIBS  = -L/usr/lib64/lib -lX11 
#X11LIBS  = -L/usr/lib/x86_64-linux-gnu/ -lX11 

CC = cc
CXX = g++

#CC=icc
#CXX=icpc


ICCFLAGS = -g -O0 -Wall -Wextra -Wno-unused-parameter
#ICCFLAGS = -g -O0 -Wall -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -Wno-write-strings
#ICCFLAGS = -g  -Wall -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -Wno-write-strings

# Set this to -pg to enable gprof profiling.
PROFILING =

FC = ifort
FFLAGS  = -g -finit-local-zero -fno-automatic 
CXXFLAGS = $(ICCFLAGS) $(UIUC_CFLAGS) $(ROOT_CFLAGS) -fPIC -Wall $(shell $$MUTRUNK/rootABI.sh $(CXX))

ALL_LIBS =  $(ROOT_LIBS) $(CERN_LIBS) $(X11LIBS)

%.o: %.C
	$(CXX) $(PROFILING) -c $(CXXFLAGS) $<
	$(CXX) -MM $(CXXFLAGS) $< > $(notdir $(patsubst %.C,%.d,$<))

%.o: %.cpp
	$(CXX) $(PROFILING) -c $(CXXFLAGS) $<
	$(CXX) -MM $(CXXFLAGS) $< > $(notdir $(patsubst %.cpp,%.d,$<))

Dictionary.cpp: MuEPair.h LinkDef.h 
	rootcint -f $@ -c $(UIUC_CFLAGS) $?

DictionaryDelayedEvent.cpp: DelayedEvent.h LinkDefDelayedEvent.h
	rootcint -f $@ -c $(UIUC_CFLAGS) $?

DictionaryLifetimeEvent.cpp: LifetimeEvent.h LifetimeEventTreeBuilder.h TLifetimeHistograms.h TLifetimeHistograms_ElectronBG.h TLifetimeHistograms_Tree.h TTPCDelayedEventHistograms.h TLifetimeHistograms_ePC.h LinkDefLifetimeEvent.h
	rootcint -f $@ -c $(UIUC_CFLAGS) $(LIFETIMEHEADERS)

DictionaryStopDistTree.cpp: stopDistPair.h stopDistPairTreeBuilder.h LinkDefStopDistTree.h 
	rootcint -f $@ -c $(UIUC_CFLAGS) $(STOPDISTHEADERS)

# Automatically include dependency information
MDINCLUDEFLAGS := $(CCINCLUDEFLAGS)
depend: $(SOURCES) $(HEADERS)
	@echo Making dependencies...
	@rm -f depend.1 depend
	@g++ -MM $(MDINCLUDEFLAGS) -D_DEPEND $(filter %.cpp,$^) > depend

libMuETree.so: $(LIB_DIR)/libMucapEvent.so MuEPair.o Dictionary.o
	$(CXX) $(PROFILING) $(ICCLIBS) --shared MuEPair.o Dictionary.o -o $@ -L$(MU) -lMucapEvent $(ROOT_LIBS)

libMuETree.rootmap: libMuETree.so LinkDef.h
ifdef ROOT6
	rootcling -f libMuETree.cxx -rml libMuETree.so -rmf libMuETree.rootmap -I$(subst :, -I,$(UIUC_INCLUDE_DIRS)) MuEPair.h ${MTA}/src/LinkDef.h
else
	rlibmap -f -o $@ -l libMuETree.so -d $(LIB_DIR)/libMucapEvent.so -c ${MTA}/src/LinkDef.h
endif

libDelayedEvent.so: DelayedEvent.o DictionaryDelayedEvent.o
	$(CXX) $(PROFILING) $(ICCLIBS) --shared DelayedEvent.o DictionaryDelayedEvent.o $(ROOT_LIBS) -o $@

libDelayedEvent.rootmap: libDelayedEvent.so LinkDefDelayedEvent.h
ifdef ROOT6
	rootcling -f libDelayedEvent.cxx -rml libDelayedEvent.so -rmf libDelayedEvent.rootmap -I$(subst :, -I,$(UIUC_INCLUDE_DIRS)) DelayedEvent.h ${MTA}/src/LinkDefDelayedEvent.h
else
	rlibmap -f -o $@ -l libDelayedEvent.so -c ${MTA}/src/LinkDefDelayedEvent.h
endif

libLifetimeEvent.so: LifetimeEvent.o DictionaryLifetimeEvent.o TLifetimeHistograms.o uiuc_tree_utils.o TLifetimeHistograms_ElectronBG.o TLifetimeHistograms_Tree.o TTPCDelayedEventHistograms.o TLifetimeHistograms_ePC.o
	$(CXX) $(PROFILING) $(ICCLIBS) --shared LifetimeEvent.o TLifetimeHistograms.o uiuc_tree_utils.o TLifetimeHistograms_ElectronBG.o TLifetimeHistograms_Tree.o TTPCDelayedEventHistograms.o TLifetimeHistograms_ePC.o DictionaryLifetimeEvent.o $(LIB_DIR)/libMucapEvent.so -o $@

libLifetimeEvent.rootmap: libLifetimeEvent.so LifetimeEvent.h LinkDefLifetimeEvent.h
ifdef ROOT6
	rootcling -f libLifetimeEvent.cxx -rml libLifetimeEvent.so -rmf libLifetimeEvent.rootmap -I$(subst :, -I,$(UIUC_INCLUDE_DIRS)) LifetimeEvent.h LifetimeEventTreeBuilder.h TLifetimeHistograms.h TLifetimeHistograms_ElectronBG.h TLifetimeHistograms_Tree.h TTPCDelayedEventHistograms.h TLifetimeHistograms_ePC.h ${MTA}/src/LinkDefLifetimeEvent.h
else
	rlibmap -f -o $@ -l libLifetimeEvent.so -d $(LIB_DIR)/libMucapEvent.so -c ${MTA}/src/LinkDefLifetimeEvent.h
endif

libStopDistTree.so: stopDistPair.o DictionaryStopDistTree.o
	$(CXX) $(PROFILING) $(ICCLIBS) --shared stopDistPair.o DictionaryStopDistTree.o $(ROOT_LIBS) -o $@

libStopDistTree.rootmap: libStopDistTree.so stopDistPair.h LinkDefStopDistTree.h
ifdef ROOT6
	rootcling -f libStopDistTree.cxx -rml libStopDistTree.so -rmf libStopDistTree.rootmap -I$(subst :, -I,$(UIUC_INCLUDE_DIRS)) stopDistPair.h stopDistPairTreeBuilder.h ${MTA}/src/LinkDefStopDistTree.h
else
	rlibmap -f -o $@ -l libStopDistTree.so -c ${MTA}/src/LinkDefStopDistTree.h
endif

mta : $(OBJECTS) $(LIB_DIR)/libMucapEvent.so $(LIB_DIR)/libTPCPulseFinder.so libMuETree.so libDelayedEvent.so libLifetimeEvent.so libStopDistTree.so
	@echo Linking...
	$(CXX) $(PROFILING) $(ICCFLAGS) -o $@ $(OBJECTS) $(LIB_DIR)/libMucapEvent.so $(LIB_DIR)/libTPCPulseFinder.so libMuETree.so libDelayedEvent.so libLifetimeEvent.so libStopDistTree.so libMinuit.so.6.24 $(ROOT_LIBS) 

clean: 
	rm -f *.d *.o Dictionary.cpp Dictionary.h DictionaryDelayedEvent.cpp DictionaryDelayedEvent.h DictionaryLifetimeEvent.cpp DictionaryLifetimeEvent.h DictionaryStopDistTree.cpp DictionaryStopDistTree.h libMuETree.so libMuETree.rootmap libDelayedEvent.so libDelayedEvent.rootmap libLifetimeEvent.so libLifetimeEvent.rootmap libStopDistTree.so libStopDistTree.rootmap mta

# Use "make cleandisplay && make" to recompile for event display use.
cleandisplay: 
	rm -f mta EventDisplaySandbox.o main.o EventDisplayTrigger.o

blank:
	touch blank

include blank $(wildcard *.d)
