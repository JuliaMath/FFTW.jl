# In MKL 2024, Intel dropped support for macOS.
# So, in CI, for the macOS jobs, we force MKL 2023 to be installed (instead of
# MKL 2024).

import TOML

using Test: @test

function main()
    root_dotci_dir = @__DIR__
    root_dir = dirname(root_dotci_dir)
    root_project_toml_filename = joinpath(root_dir, "Project.toml")
    project = TOML.parsefile(root_project_toml_filename)
    old_mkl_jll_compat = project["compat"]["MKL_jll"]
    old_mkl_jll_compat_list = strip.(split(strip(old_mkl_jll_compat), ","))

    # Regression test to catch if we ever drop support for MKL_jll 2023.
    # Because, if we ever do choose to drop support for MKL_jll 2023, then this
    # entire script should probably be deleted, since Intel dropped support for
    # macOS in MKL 2024.
    @test "2023" in old_mkl_jll_compat_list

    new_compat = "2023"
    project["compat"]["MKL_jll"] = new_compat # force MKL_jll 2023
    open(root_project_toml_filename, "w") do io
        TOML.print(io, project)
    end
    @info "Changed the compat entry for MKL_jll to: $(new_compat)"
    return nothing
end

main()
