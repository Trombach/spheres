CC = clang++-6.0
CCFLAGS = -O3 -std=c++11 -fopenmp=libiomp5 
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
MATLIB = ~/dlib-19.2/

SRCS := $(wildcard *.cc)
MAIN := optimize.cc analyze.cc match.cc getxyz.cc ico-subgraph.cc global.cc eigenvectorFollowing.cc $(wildcard app-*)
APPS := $(wildcard apps/*.cc)
EXES := $(patsubst %.cc, %, $(MAIN))
EXES += $(patsubst %.cc, %, $(APPS))
SRCS := $(filter-out $(MAIN), $(SRCS))

LFLAGS = -lconfig++ -llapack
INCLUDE =-I $(MATLIB)
SRCS := $(filter-out gsl%, $(SRCS))

OBJS = $(patsubst %.cc, %.o, $(SRCS))
$(info SRCS = $(SRCS))
$(info OBJS = $(OBJS))
$(info EXES = $(EXES))
$(info APPS = $(APPS))

OBJDIR = objects
DEPDIR = .dep
$(shell mkdir -p $(OBJDIR) >/dev/null)
$(shell mkdir -p $(DEPDIR) >/dev/null)
$(shell mkdir -p $(DEPDIR)/apps >/dev/null)

FLAGS = $(CCFLAGS) $(LFLAGS) $(INCLUDE)
COMPILE.cc = $(CC) $(DEPFLAGS) $(CCFLAGS) $(INCLUDE)
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d


.PHONY: all
all: $(EXES)

$(EXES): %: $(OBJDIR)/%.o $(addprefix $(OBJDIR)/,$(OBJS))
	$(CC) $(FLAGS) -o $@ $^

%.o: %.cc
$(OBJDIR)/%.o: %.cc
$(OBJDIR)/%.o: %.cc | $(DEPDIR)/%.d
	$(COMPILE.cc) -c -o $@ $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

.PHONY: distclean
distclean:
	rm $(EXES) $(OBJDIR)/*.o $(DEPDIR)/*.d $(OBJDIR)/apps/*.o $(DEPDIR)/apps/*.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))
-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(MAIN)))
-include $(patsubst %,$(DEPDIR)/apps/%.d,$(basename $(APPS)))
