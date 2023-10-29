.PHONY: all, update, instantiate, test, clean, pages, format-check, format

# Julia command
JULIA = julia
PKG_NAME = ONSAS

# Update dependecies
update:
	$(JULIA) --project=. -e 'import Pkg; Pkg.update()'  && $(JULIA) --project=./docs/ -e 'import Pkg; Pkg.update()'

# Update dependecies
instantiate:
	$(JULIA) --project=. -e 'import Pkg; Pkg.instantiate()'

# Run all tests
test:
	$(JULIA) --project=. -e 'include("test/runtests.jl")'

# Make docs
pages: 
	$(JULIA) --project=./docs/ -e 'import Pkg; Pkg.instantiate(); include("./docs/make.jl")'

# Remove generated files
clean:
	find . -name "*.msh" -type f -exec rm -f {} +
	find . -name "*.vti" -type f -exec rm -f {} +
	find . -name "*.vtu" -type f -exec rm -f {} +

# Check JuliaFormatter is applied
format-check:
	$(JULIA) --project=. -e '\
	using JuliaFormatter; \
	files = read(`find . -name "*.jl"`, String); \
	for file in split(files, "\n"); \
		if file != "" && !isdir(file) && !format(file, verbose=true, overwrite=false); \
			println("Formatting issue found in: ", file); \
			exit(1); \
		end; \
	end; \
	println("No Formatting issues found.")'

# Apply JuliaFormatter is applied
format:
	$(JULIA) --project=. -e 'using JuliaFormatter; format(".", overwrite=true, verbose=true)'
