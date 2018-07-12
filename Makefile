# Stores root folder
TOP_DIR=$(shell git rev-parse --show-toplevel)

# Colours
RED = "\\033[31m"
YELLOW = "\\033[33m"
GREEN = "\\033[92m"
NC = "\\033[0m"
RV = "\\033[7m"
SB = "\\033[5m"

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
	@printf "$(GREEN)*** Generating documentation ***$(NC)\n"
	@printf "\n"
	Rscript -e 'devtools::document()'
	
## install : install the package
.PHONY: install
install: document
	@printf "\n"
	@printf "$(GREEN)*** Installing package ***$(NC)\n"
	@printf "\n"
	Rscript -e 'devtools::install()'
	@printf "\n"

## clean : remove man folder and test output images
.PHONY: clean
clean:
	@printf "\n"
	@printf "$(YELLOW)*** Cleaning ***$(NC)\n"
	@printf "\n"
	rm -fr man/
	rm -f sandbox/Rplots.pdf
	@printf "\n"
	
## uninstall : uninstall the package
.PHONY: uninstall
uninstall:
	@printf "\n"
	@printf "$(RED)*** Uninstalling package ***$(NC)\n"
	@printf "\n"
	Rscript -e "tryCatch(remove.packages('whereipm'), error = function(e) {})"
	@printf "\n"

## test : run test on package
.PHONY: test
test:
	@printf "\n"
	@printf "$(GREEN)$(FR)*** Running test ***$(NC)\n"
	@printf "\n"
	Rscript -e "devtools::test()"
	@printf "\n"
	
## coverage : determine test coverage
.PHONY: coverage
coverage:
	@printf "\n"
	@printf "$(GREEN)*** Determining test coverage ***$(NC)\n"
	@printf "\n"
	@mkdir -p docs/coverage/
	@Rscript -e "library(covr); cov = package_coverage(); report(cov, \"$(TOP_DIR)/docs/coverage/coverage.html\"); print(paste(\"Coverage_percent: --\", round(percent_coverage(cov), 2), \"--\"))"
	
## rehearse : unistall package, clean, install again and run test
.PHONY: rehearse
rehearse: uninstall clean install test
	
