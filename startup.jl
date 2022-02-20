# activate env based on CWD
using Pkg
proj = Base.current_project()
proj === nothing || Pkg.activate(proj)
# allows for running project specific startup code in the Project.toml entry [juliaenv] startup = "..."
using ProjectX
