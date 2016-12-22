all:
	f2py -c --compiler=unix --f90flags="-O3 -fopenmp -ffixed-line-length-0" -m UDP_modules critfile.f90 UDP_modules.f90	
#	f2py -c --debug-capi --compiler=unix --f90flags="-fbounds-check -Wall" -m UDP_modules critfile.f90 UDP_modules.f90	

clean:
	rm -f UDP_modules.so UDPClust.pyc
