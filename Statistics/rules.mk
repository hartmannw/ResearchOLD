# Rules for Statistics directory.
local_dir  := Statistics
local_relsrc  := DiagonalGaussian.cc HiddenMarkovModel.cc HmmSet.cc \
	MixtureOfDiagonalGaussians.cc PosteriorgramGenerator.cc
local_src  := $(addprefix $(local_dir)/,$(local_relsrc))
local_relexec  := test_posteriorgram
local_exec  := $(addprefix $(local_dir)/,$(local_relexec))

SRCS += $(local_src)
EXECS += $(local_exec)
