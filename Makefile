.PHONY: install run test
install:
	Rscript install.R
run:
	Rscript run.R
test:
	Rscript tests/run.R
