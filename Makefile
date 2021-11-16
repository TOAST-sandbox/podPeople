install: ngs-tools/tools/tax/Makefile
	cd ngs-tools/tools/tax && ./quickbuild.sh
	@echo "remember: export PATH=$$(pwd -P)/ngs-tools/tools/tax/bin"

ngs-tools/tools/tax/Makefile:
	git clone https://github.com/ncbi/ngs-tools.git --branch tax


