iteraciones := 20000
quemados := 5000

TengoSuenho.pdf : datos.dat graficando.py
	python graficando.py
	evince TengoSuenho.pdf &

datos.dat : mcmc_lotkavolterra.x
	./mcmc_lotkavolterra.x 20000 100 > datos.dat

mcmc_lotkavolterra.x : diferencial.c
	cc diferencial.c -lgsl -lgslcblas -lm -o mcmc_lotkavolterra.x

clean:
	rm -f mcmc_lotkavolterra.x datos.dat TengoSuenho.pdf