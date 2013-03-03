tarname:=tree19
ARCH=x86_64

# Make.$(ARCH) sets many of the variables that are then used in
# Make.generic.  Leaving it out can cause problems, for example,
# Make.$(ARCH) wants to redefine 'TAR'.  On the other hand, leaving
# it in can cause problems (apparently) with older versions of GNUmake
# because, e.g., assignements like LOADLIBES:=-foo $(LOADLIBES) get
# 'executed' multiple times.  Once in this Makefile, and once in the
# daughter makefiles.

###include Make-common/Make.$(ARCH)

# we should also go into SDFcvt, SDF2fld, commtst, lsv, lsvtst, (anything
# else?) and build them.  Possibly under a different target?

all: All

# Make.generic has targets for clean, all, etc., so we need to
# spell them a little differently in this file...
include Make-common/Make.generic

subdirs:= libSDF libmpmy libsw apps/SDFcvt

All:
	for dir in $(subdirs); do (cd $$dir; $(MAKE) ARCH=$(ARCH) all); done

Depends :
	for dir in $(subdirs); do (cd $$dir; $(MAKE) ARCH=$(ARCH) depends); done

Clean : 
	for dir in $(subdirs); do (cd $$dir; $(MAKE) ARCH=$(ARCH) clean); done

# $(treedir)
