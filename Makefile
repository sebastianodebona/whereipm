# Colours
GREEN = "\\033[92m"
NC = "\\033[0m"

# help
.PHONY: help
help: Makefile
	@printf "\n"
	@printf "$(GREEN)*** Please use 'make <target>' where <target> is one of: ***$(NC)\n"
	@printf "\n"
	@sed -n 's/^## /    /p' $< | column -t -s ":"
	@printf "\n"

## document : generate documentation with roxygen
.PHONY: document
document:
	@printf "\n"
	@printf "Generating documentation\n"
	@printf "\n"
	Rscript -e 'devtools::document()'
	
## install : install the package
.PHONY: install
install: document
	Rscript -e 'devtools::install()'
	
## clean : remove man folder
.PHONY:
clean:
	rm -fr man/
	