"""
    create_run(run_name::AbstractString)

Generates a generic input file ("run_name.in") loaded with options for the run in the current directory.
"""
function create_run(run_name::AbstractString)
    open(string(run_name,".in"),"w") do io
    println(io,
        "run_type 1\nsites 1:36\nsite_list sitelist.txt\nAO 1\nradius 2.2\ncontinue_from lastrun.mat")
    end
end