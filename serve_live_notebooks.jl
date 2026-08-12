#!/usr/bin/env julia

include(joinpath(@__DIR__, "deployment_notebooks.jl"))
using .DeploymentNotebooks
import PlutoSliderServer

function usage()
    println("Usage: julia --project=pluto-deployment-environment serve_live_notebooks.jl [--host HOST] [--port PORT]")
end

function options(args)
    host = "127.0.0.1"
    port = 1234
    i = 1
    while i <= length(args)
        if args[i] == "--host" && i < length(args)
            host = args[i + 1]
            i += 2
        elseif args[i] == "--port" && i < length(args)
            port = parse(Int, args[i + 1])
            i += 2
        elseif args[i] in ("-h", "--help")
            usage()
            exit()
        else
            usage()
            error("Unknown or incomplete option: $(args[i])")
        end
    end
    (; host, port)
end

settings = options(ARGS)
notebooks = notebook_paths()
isempty(notebooks.live) && error("live-notebooks.yml does not contain a live notebook")

println("Starting live server for: $(join(notebooks.live, ", "))")
PlutoSliderServer.run_directory(
    SOURCE_DIRECTORY;
    notebook_paths=notebooks.live,
    SliderServer_host=settings.host,
    SliderServer_port=settings.port,
    SliderServer_watch_dir=false,
    SliderServer_serve_static_export_folder=false,
    Export_enabled=false,
)
