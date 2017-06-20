using BinDeps
using BinDeps: builddir

BinDeps.@setup

const FFTW_VER = v"3.3.6-pl1"

if is_windows()
    const libfftw_name = "libfftw3-3"
    const libfftwf_name = "libfftw3f-3"
else
    const libfftw_name = "libfftw3_threads"
    const libfftwf_name = "libfftw3f_threads"
end

libfftw = library_dependency(libfftw_name, aliases=[replace(libfftw_name, "3", "")])
libfftwf = library_dependency(libfftwf_name, aliases=[replace(libfftwf_name, "3", "")])

provides(AptGet, "libfftw3-double3", [libfftw], os=:Linux)
provides(AptGet, "libfftw3-single3", [libfftwf], os=:Linux)
provides(Pacman, "fftw", [libfftw, libfftwf], os=:Linux)
provides(Zypper, "libfftw3", [libfftw, libfftwf], os=:Linux)
provides(Yum, "fftw", [libfftw, libfftwf], os=:Linux)
provides(BSDPkg, "fftw3", [libfftw, libfftwf], os=:FreeBSD)

if is_windows()
    using WinRPM
    provides(WinRPM.RPM, "fftw3", [libfftw, libfftwf], os=:Windows)
elseif is_apple()
    using Homebrew
    provides(Homebrew.HB, "fftw", [libfftw, libfftwf], os=:Darwin)
end

general_config = ["--prefix=" * abspath(builddir(libfftw)),
                  "--libdir=" * abspath(libdir(libfftw)),
                  "--bindir=" * abspath(bindir(libfftw))]

fftw_config = ["--enable-shared", "--disable-fortran", "--disable-mpi", "--enable-threads"]
fftw_enable_single = "--enable-single"

if Sys.ARCH === :ppc
    append!(fftw_config, ["--enable-altivec", "--enable-fma"])
elseif Sys.ARCH === :x86_64
    append!(fftw_config, ["--enable-sse2", "--enable-fma"])
end

if is_windows()
    append!(fftw_config, ["--with-our-malloc", "--with-combined-threads"])
    Sys.ARCH === :x86_64 || push!(fftw_config, "--with-incoming-stack-boundary=2")
end

provides(Sources, URI("http://www.fftw.org/fftw-$FFTW_VER.tar.gz"), [libfftw, libfftwf])

provides(BuildProcess, (@build_steps begin
    GetSources(libfftw)
    CreateDirectory(builddir(libfftw))
    @build_steps begin
        ChangeDirectory(builddir(libfftw))
        FileRule(joinpath(libdir(libfftw), libfftw_name * "." * Libdl.dlext), @build_steps begin
            CreateDirectory(libdir(libfftw))
            `$(joinpath(srcdir(libfftw), "fftw-$FFTW_VER", "configure")) $general_config $fftw_config`
            `$MAKE_CMD`
            `$MAKE_CMD install`
        end)
        FileRule(joinpath(libdir(libfftw), libfftwf_name * "." * Libdl.dlext), @build_steps begin
            `$(joinpath(srcdir(libfftw), "fftw-$FFTW_VER", "configure")) $general_config $fftw_config $fftw_enable_single`
            `$MAKE_CMD`
            `$MAKE_CMD install`
        end)
    end
end), [libfftw, libfftwf])

if is_windows()
    BinDeps.@install Dict([:libfftw3 => :libfftw, :libfftw3f => :libfftwf])
else
    BinDeps.@install Dict([:libfftw3_threads => :libfftw, :libfftw3f_threads => :libfftwf])
end
