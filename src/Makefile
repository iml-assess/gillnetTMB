all: gillnetTMB.so

armatrix.so: gillnetTMB.cpp
	Rscript --vanilla -e "TMB::compile('gillnetTMB.cpp','-O0 -g')"

clean:
	rm -rf *o