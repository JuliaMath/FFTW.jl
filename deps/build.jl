using BinDeps
using BinDeps: builddir

BinDeps.@setup

const FFTW_VER = v"3.3.6-pl1"

if is_windows()
    const libfftw_name = "libfftw3"
    const libfftwf_name = "libfftw3f"
else
    const libfftw_name = "libfftw3_threads"
    const libfftwf_name = "libfftw3f_threads"
end

libfftw = library_dependency(libfftw_name, aliases=[replace(libfftw_name, "3", "")])
libfftwf = library_dependency(libfftwf_name, aliases=[replace(libfftwf_name, "3", "")])

provides(AptGet, "libfftw3", [libfftw, libfftwf], os=:Linux)
provides(Pacman, "fftw", [libfftw, libfftwf], os=:Linux)
provides(Zypper, "libfftw3", [libfftw, libfftwf], os=:Linux)
provides(Yum, "fftw", [libfftw, libfftwf], os=:Linux)
provides(BSDPkg, "fftw3", [libfftw, libfftwf], os=:FreeBSD)

if is_windows()
    using WinRPM
    provides(WinRPM.RPM, "fftw", [libfftw, libfftwf], os=:Windows)
elseif is_apple()
    using Homebrew
    provides(Homebrew.HB, "fftw", [libfftw, libfftwf], os=:Darwin)
end

general_config = ["--prefix=" * abspath(builddir(libfftw)),
                  "--libdir=" * abspath(libdir(libfftw)),
                  "--bindir=" * abspath(bindir(libfftw))]

try
    # There might not be a cc defined or available in the path
    machine = readchomp(`cc -dumpmachine`)
    push!(general_config, "--build=" * machine)
end

fftw_config = ["--enabled-shared", "--disable-fortran", "--disable-mpi", "--enable-threads"]
fftw_enable_single = "--enable-single"

if Sys.ARCH === :ppc
    push!(fftw_config, "--enable-altivec")
elseif Sys.ARCH === :x86_64
    append!(fftw_config, ["--enable-sse2", "--enable-fma"])
end

if is_windows()
    append!(fftw_config, ["--with-our-malloc", "--with-combined-threads"])
    Sys.ARCH === :x86_64 || push!(fftw_config, "--with-incoming-stack-boundary=2")
end

const MAKE = is_bsd() && !is_apple() ? `gmake` : `make`

provides(Sources, URI("http://www.fftw.org/fftw-$FFTW_VER.tar.gz"), [libfftw, libfftwf])

provides(BuildProcess, (@build_steps begin
    GetSources(libfftw)
    CreateDirectory(joinpath(builddir(libfftw), "libfftw"))
    @build_steps begin
        ChangeDirectory(joinpath(builddir(libfftw), "libfftw"))
        FileRule(joinpath(libdir(libfftw), "libfftw." * Libdl.dlext), @build_steps begin
            CreateDirectory(libdir(libfftw))
            `$(joinpath(".", "configure")) $general_config $fftw_config $fftw_enable_single`
            `$MAKE`
        end)
    end
end), libfftw)

BinDeps.@install Dict([:libfftw => :libfftw, :libfftwf => :libfftwf])
