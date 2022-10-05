mqcroot = "$(HOME)/.mqc"
FC = gfortran
ifeq ($(FC),gfortran)
        FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
        MQCOBJS = -I$(mqcroot)/GNU/mod
        LIBS = -llapack -lblas $(mqcroot)/GNU/lib/libmqc.a
else ifeq ($(FC),pgfortran)
        FCFLAGS = -Mallocatable=03 -r8 -i8 -mp
        MQCOBJS = -module $(mqcroot)/PGI/mod
        LIBS = -llapack -lblas $(mqcroot)/PGI/lib/libmqc.a
endif

all: compile

compile: component.f03
	$(FC) $(FCFLAGS) $(MQCOBJS) -o component.exe component.f03 $(LIBS) 

clean:
	rm -f component.exe component.o 
