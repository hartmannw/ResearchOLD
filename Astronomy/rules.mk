# Specific make rules for the Astronomy directory
local_dir  := Astronomy
local_relsrc  := OrbitalBody.cc GravitationalSystem.cc
local_src  := $(addprefix $(local_dir)/,$(local_relsrc))
local_relexec  := ModelSolarSystem test
local_exec  := $(addprefix $(local_dir)/,$(local_relexec))

SRCS += $(local_src)
EXECS += $(local_exec)
