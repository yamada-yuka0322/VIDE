# Make.$(ARCH) sets many of the variables used below including:
# CC, CFLAGS, AS, RANLIB, objdir, objsuf, asmdir

treedir=..
treedir_sed=\.\.
appexcludes=
# Put everything into the libsw library!!!
libname=libsw

# The "SDF" library
src:= \
SDF-parse.c 	SDFfuncs.c	SDFget.c	SDFhdrio.c

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

# These rules are slight modifications of the builtin ones
%.c : %.y
	-($(YACC.y) $<; mv -f y.tab.c $@)

%.c : %.l
#  commands to execute (built-in):
	-($(LEX.L) $< > y.lex.c; $(RM) $@; mv y.lex.c $@)

# Makedepends can't pick these up...
$(objdir)/SDF-parse$(objsuf): SDF-private.h SDF-lex.c SDF-parse.c
	(cd $(objdir); $(CC) $(CFLAGS) -I../../../include/libsdf -I../.. -c ../../SDF-parse.c)

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf): $(treedir)/include/libsdf/protos.h $(treedir)/include/libsdf/Msgs.h
$(objdir)/SDF-parse$(objsuf): $(treedir)/include/libsdf/gccextensions.h SDF-private.h stdio.h
$(objdir)/SDF-parse$(objsuf): $(treedir)/include/libsdf/mpmy_io.h
$(objdir)/SDF-parse$(objsuf): $(treedir)/include/libsdf/SDF.h $(treedir)/include/libsdf/Malloc.h
$(objdir)/SDF-parse$(objsuf): $(treedir)/include/libsdf/error.h SDF-lex.c
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDF-parse$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFfuncs$(objsuf): $(treedir)/include/libsdf/Msgs.h $(treedir)/include/libsdf/gccextensions.h
$(objdir)/SDFfuncs$(objsuf): stdio.h $(treedir)/include/libsdf/mpmy_io.h SDF-private.h
$(objdir)/SDFfuncs$(objsuf): $(treedir)/include/libsdf/SDF.h
$(objdir)/SDFfuncs$(objsuf): $(treedir)/include/libsdf/byteswap.h $(treedir)/include/libsdf/Malloc.h
$(objdir)/SDFfuncs$(objsuf): $(treedir)/include/libsdf/error.h $(treedir)/include/libsdf/protos.h
$(objdir)/SDFfuncs$(objsuf): $(treedir)/include/libsdf/mpmy.h $(treedir)/include/libsdf/timers.h
$(objdir)/SDFfuncs$(objsuf):
$(objdir)/SDFget$(objsuf): stdio.h $(treedir)/include/libsdf/mpmy_io.h
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf):
$(objdir)/SDFget$(objsuf): $(treedir)/include/libsdf/error.h
$(objdir)/SDFget$(objsuf): $(treedir)/include/libsdf/gccextensions.h $(treedir)/include/libsdf/SDF.h
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf): $(treedir)/include/libsdf/Msgs.h $(treedir)/include/libsdf/gccextensions.h
$(objdir)/SDFhdrio$(objsuf): $(treedir)/include/libsdf/error.h SDF-private.h
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf): stdio.h $(treedir)/include/libsdf/mpmy_io.h
$(objdir)/SDFhdrio$(objsuf):
$(objdir)/SDFhdrio$(objsuf): $(treedir)/include/libsdf/SDF.h
