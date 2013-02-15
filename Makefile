# A non-recursive makefile
# Initialize variables used by submakefiles.
SRCS :=
EXECS :=

OBJS = $(subst .cc,.o,$(SRCS))
DEPS = $(subst .cc,.d,$(SRCS))

# Standard flags and programs
CC=g++
CPPFLAGS=-Wall -pedantic -O2 -std=c++11
DEPFLAGS=-MD -MP -MF
LDFLAGS=
LDLIBS=
RM := rm -f

# Handle directory structure
INCLUDE_DIRS := Utilities Statistics FileUtilities Astronomy \
	AcousticUnitDiscovery
CPPFLAGS     += $(addprefix -I ,$(INCLUDE_DIRS))
vpath %.h $(INCLUDE_DIRS)

# Initialize the all target here so that it is the default.
all:

# Include all sub makefiles
include Utilities/rules.mk
include Statistics/rules.mk
include FileUtilities/rules.mk
include Astronomy/rules.mk
include AcousticUnitDiscovery/rules.mk

.PHONY: all
all: $(EXECS)

.PHONY: clean
clean:
	$(RM) $(OBJS) $(EXECS) $(DEPS)


# Only include dependencies if we are not cleaning repository.
ifneq "$(MAKECMDGOALS)" "clean"
  include $(dependencies)
endif

# link
${EXECS} : $(OBJS)
	$(CC) -o $@ $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $(OBJS) $@.cc

# compile .o and generate dependencies
%.o : %.cc
	$(CC) -c -o $@ $(CPPFLAGS) $(DEPFLAGS) $(@:.o=.d) $<
