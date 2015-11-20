iteraciones := 200000
quemados := 50000

graficas.py: datos.dat
	python graficas.py

datos.dat: mcmc_solar.x
	./mcmc_solar.x $(iteraciones) $(quemados) > datos.dat 



mcmc_solar.x: retro.c
	cc retro.c -lgsl -lgslcblas -lm -o mcmc_solar.x
clean:
	rm -f datos.dat mcmc_solar.x