# Specific make rules for the AcousticUnitDiscovery directory
local_dir  := AcousticUnitDiscovery
local_relsrc  := DynamicTimeWarp.cc MultiBestPath.cc
local_src  := $(addprefix $(local_dir)/,$(local_relsrc))
local_relexec  := CreateSimilarityMatrix GeneratePronunciations testdtw \
	testmultibest
local_exec  := $(addprefix $(local_dir)/,$(local_relexec))

SRCS += $(local_src)
EXECS += $(local_exec)
