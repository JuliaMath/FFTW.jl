using BinDeps
using BinDeps: builddir, depsdir, libdir

# Binaries is not a recognized provider on Linux >:/
modified_defaults = false
if !in(BinDeps.Binaries, BinDeps.defaults)
    unshift!(BinDeps.defaults, BinDeps.Binaries)
    modified_defaults = true
end

BinDeps.@setup

const FFTW_VER = v"3.3.6-pl2"

if Sys.iswindows()
    const libfftw_name = "libfftw3"
    const libfftwf_name = "libfftw3f"
else
    const libfftw_name = "libfftw3_threads"
    const libfftwf_name = "libfftw3f_threads"
end

# Why can't everyone just agree on what to call this library...
function makealiases(lib)
    major = string(FFTW_VER.major)
    nover = replace(lib, major, "")
    return String[
        nover,
        join([lib, Libdl.dlext, major], "."),
        join([nover, Libdl.dlext, major], "."),
        lib * "-" * major,
        nover * "-" * major,
    ]
end

libfftw = library_dependency(libfftw_name, aliases=makealiases(libfftw_name))
libfftwf = library_dependency(libfftwf_name, aliases=makealiases(libfftwf_name))

const URL = "https://github.com/ararslan/fftw-builder/releases/download/v$FFTW_VER/libfftw-$FFTW_VER"

# Mapping of Sys.MACHINE to (url, sha) for precompiled binaries from fftw-builder
const downloads = Dict(
    "x86_64-pc-linux-gnu" => ("$URL-linux-x86_64.tar.gz",
                              "f8382542f93b391590faee540d71ed8022d5e04037ae664e8f8a6a445b2ea876"),
    "i686-pc-linux-gnu"   => ("$URL-linux-i686.tar.gz",
                              "b3a7e6726470de92407259ad198d9d9e256df77a3e0276de5d9a81965f5c5176"),
    "x86_64-apple-darwin" => ("$URL-osx-x86_64.tar.gz",
                              "97fa4cca555d4587c1d2a02d749c92c513cc59eb209a19d3882009006d911c6b"),
    "x86_64-w64-mingw32"  => ("$URL-win-x86_64.zip",
                              "40199c4de73e23c0009acb625b14999130d489621d357c9d88003559089d69eb"),
    "i686-w64-mingw32"    => ("$URL-win-i686.zip",
                              "47eeac9c8a8a869ac0d9114e4a57dc6463aa4dd4e7dd12eaf73aea9993c04b8f"),
)

const machine = Sys.isapple() ? "x86_64-apple-darwin" : Sys.MACHINE

if haskey(downloads, machine)
    url, sha = downloads[machine]
    provides(Binaries, URI(url), [libfftw, libfftwf], SHA=sha, os=BinDeps.OSNAME,
             unpacked_dir=joinpath("usr", "lib"), installed_libpath=libdir(libfftw))
    scratch = false
elseif Sys.KERNEL === :FreeBSD
    provides(BSDPkg, "fftw3", [libfftw, libfftwf], os=:FreeBSD)
    scratch = false
else
    info("No precompiled binaries found for your system. Building from scratch...")
    scratch = true
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

if Sys.iswindows()
    append!(fftw_config, ["--with-our-malloc", "--with-combined-threads"])
    Sys.ARCH === :x86_64 || push!(fftw_config, "--with-incoming-stack-boundary=2")
end

# Make it harder to build from scratch
if scratch
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
end

if Sys.iswindows()
    BinDeps.@install Dict([:libfftw3 => :libfftw, :libfftw3f => :libfftwf])
else
    BinDeps.@install Dict([:libfftw3_threads => :libfftw, :libfftw3f_threads => :libfftwf])
end

if modified_defaults
    shift!(BinDeps.defaults)
end
