# Specific make rules for the Utilities directory
local_dir  := Utilities
local_relsrc  := StringFunctions.cc Matrix.cc MatrixFunctions.cc
local_src  := $(addprefix $(local_dir)/,$(local_relsrc))
local_relexec  := testmatrix
local_exec  := $(addprefix $(local_dir)/,$(local_relexec))

SRCS += $(local_src)
EXECS += $(local_exec)
