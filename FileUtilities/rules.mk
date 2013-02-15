# Specific make rules for the Utilities directory
local_dir  := FileUtilities
local_relsrc  := SpeechFeatures.cc ImageIO.cc
local_src  := $(addprefix $(local_dir)/,$(local_relsrc))
local_relexec  := ConvertCep2Ascii ConvertCep2Htk test
local_exec  := $(addprefix $(local_dir)/,$(local_relexec))

SRCS += $(local_src)
EXECS += $(local_exec)
