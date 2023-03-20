using Pkg
Pkg.activate(Base.current_project())
Pkg.instantiate()
using Conda
Conda.add("POT")
ENV["PYTHON"] = ""
Pkg.build("PyCall")