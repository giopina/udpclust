all:
	f2py -c -m UDP_modules critfile.f90 UDP_modules.f90	

clean:
	rm -f UDP_modules.so
